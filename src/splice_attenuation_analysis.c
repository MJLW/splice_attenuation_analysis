#include <cpliceai/range.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <cpliceai/predict.h>
#include <cpliceai/utils.h>

#include "bcftools/gff.h"
#include "bcftools/regidx.h"

#include "logging/log.h"

#define TENSORFLOW_LOG_LEVEL_SILENT 1 // Only warnings and above
#define SPLICE_SITE_MALLOC_START_COUNT 1000000

#define SPLICEAI_CONTEXT_PADDING 5000
#define SPLICEAI_WINDOW_PADDING 50

#define BASES "ACGT"

#define COMPACT_HEADER_STRING "chr\tpos\tgene\ttid\tstrand\ttype\tvar_pos\tvar_ref\tvar_alt\tvar_dp\tcanonical_ref\tcanonical_alt\tpair_ref\tpair_alt\tcompetitor_ref\tcompetitor_alt\n"
#define PREFIX_HEADER_STRING "chr\tpos\tgene\ttid\tstrand\ttype"

typedef enum {
    ACCEPTOR = ACCEPTOR_POS,
    DONOR = DONOR_POS
} SpliceSiteType;

typedef struct {
    const char *chr;
    uint32_t rid; // Chromosome/region encoding, for quick comparisons
    uint64_t pos;
    SpliceSiteType type;

    const char *gene_name;
    uint64_t gene_beg;
    uint64_t gene_end;
    uint32_t tid;
    int strand;

    uint64_t has_pair:1, pair_position:63;
    uint64_t has_competitor:1, competitor_position:63;

    // uint32_t predicted;
    // float score[NUM_SCORES]; // Should only be used if predicted==TRUE
} SpliceSite;

typedef struct {
    uint32_t n, m;
    SpliceSite *a;
} SpliceSites;

SpliceSite init_splice_site(const char *chr, const uint32_t rid, const uint64_t pos, const SpliceSiteType type, const char *gene_name, const uint64_t beg, const uint64_t end, const uint32_t tid, const int strand) {
    return (SpliceSite) { chr, rid, pos, type, gene_name, beg, end, tid, strand, 0, 0, 0, 0};
}

typedef struct {
    const char *chr;
    uint64_t pos;
    const char *gene_name;
    uint32_t tid;
    int strand;
    SpliceSiteType type;
    uint64_t var_pos;
    float var_ref;
    float var_alt;
    int var_dp;
    float canonical_ref, canonical_alt;
    int has_pair;
    float pair_ref, pair_alt;
    int has_competitor;
    float competitor_ref, competitor_alt;
} CompactOutputRow;

typedef enum {
    REF,
    ALT
} ProbabilityType;

typedef struct {
    const char *chr;
    uint64_t pos;
    const char *gene_name;
    uint32_t tid;
    int strand;
    SpliceSiteType type;
    ProbabilityType prob_type;
    float *probabilities;
} ProbabilityOutputRow; 

void splice_sites_init(const uint32_t m, SpliceSites *sites) {
    sites->n = 0;
    sites->m = m;
    sites->a = malloc(m * sizeof(SpliceSite));
}

void splice_sites_destroy(SpliceSites *sites) {
    free(sites->a);
    sites->a = NULL;
    sites->n = 0;
    sites->m = 0;
}

int find_next_site_with_type(const SpliceSites sites, const int tid, const int start, const int end, const int direction, const SpliceSiteType target_type, SpliceSite *out) {
    if (direction > 0) {
        for (int i = start+1; i < end; i++) {
            SpliceSite site = sites.a[i];
            if (site.tid != tid) break;
            if (site.type != target_type) continue;

            *out = site;
            return 0;
        }
    } else {
        for (int i = start-1; i >= end; i--) {
            SpliceSite site = sites.a[i];
            if (site.tid != tid) break;
            if (site.type != target_type) continue;

            *out = site;
            return 0;
        }
    }

    return 1;
}

SpliceSites get_splice_sites_from_gff(const gff_t *gff) {
    SpliceSites sites;
    splice_sites_init(SPLICE_SITE_MALLOC_START_COUNT, &sites);

    regidx_t *transcripts = gff_get((gff_t *) gff, idx_tscript); // Removes const qualification for this call as it's not typed const, but it is a simple getter without consequences
    regitr_t *itr = regitr_init(transcripts);
    while (regitr_loop(itr)) {
        const gf_tscript_t *tr = regitr_payload(itr, gf_tscript_t *);

        // Only want protein coding and stranded transcripts
        if ((tr->type != GF_PROTEIN_CODING) | (tr->gene->type != GF_PROTEIN_CODING)) continue;
        if (tr->strand == STRAND_UNK) log_error("Transcript has an UNKNOWN strand. Skipping...", itr->seq);

        const char *chr = itr->seq;
        const int rid = tr->gene->iseq;

        // Only want things located on chromosomes
        if (strncmp(chr, "chr", 3) != 0) continue; // WARN: This is based on RUMC GFF files, they may not always be prefixed with chr
        if (sites.n + 2 >= sites.m) sites.a = realloc(sites.a, (sites.m *= 2) * sizeof(SpliceSite));

        // Extract stuff from pointers, no practical benefit, just looks a bit cleaner
        const char *gene_name = tr->gene->name;
        const uint64_t tr_beg = tr->gene->beg, tr_end = tr->gene->end;
        const uint32_t tid = tr->id;
        const int strand = tr->strand;

        for (int i = 0; i < tr->ncds; i++) {
            const gf_cds_t *cds = tr->cds[i];
            const uint64_t cds_beg = cds->beg;
            const uint64_t cds_end = cds->beg + cds->len - 1; // Offset by -1 to get closed end coordinate

            if (tr->strand == STRAND_FWD) {
                if (i != 0) sites.a[sites.n++] = init_splice_site(chr, rid, cds_beg, ACCEPTOR, gene_name, tr_beg, tr_end, tid, strand); // As long as its not the first exon
                if (i != tr->ncds-1) sites.a[sites.n++] = init_splice_site(chr, rid, cds_end, DONOR, gene_name, tr_beg, tr_end, tid, strand); // As long as its not the last exon
            } else if (tr->strand == STRAND_REV) {
                if (i != 0) sites.a[sites.n++] = init_splice_site(chr, rid, cds_beg, DONOR, gene_name, tr_beg, tr_end, tid, strand); // As long as its not the last exon (approaches from 3' -> 5')
                if (i != tr->ncds-1) sites.a[sites.n++] = init_splice_site(chr, rid, cds_end, ACCEPTOR, gene_name, tr_beg, tr_end, tid, strand); // As long as its not the first exon (approaches from 3' -> 5')
            }
        }

    }

    // Get pairs
    const int pair_start = 0;
    const int pair_end = sites.n;
    for (int i = 0; i < sites.n; i++) {
        SpliceSite *site = &sites.a[i];

        SpliceSite pair;
        int pair_missing;
        if (site->type == DONOR) {
            if (site->strand == STRAND_FWD) {
                pair_missing = find_next_site_with_type(sites, site->tid, i, pair_end, 1, ACCEPTOR, &pair);
            } else { // site.strand == STRAND_REV
                pair_missing = find_next_site_with_type(sites, site->tid, i, pair_start, -1, ACCEPTOR, &pair);
            }
        } else { // site.type == ACCEPTOR
            if (site->strand == STRAND_FWD) {
                pair_missing = find_next_site_with_type(sites, site->tid, i, pair_start, -1, DONOR, &pair);
            } else { // site.strand == STRAND_REV
                pair_missing = find_next_site_with_type(sites, site->tid, i, pair_end, 1, DONOR, &pair);
            }
        }

        if (pair_missing) {
            log_error("Unexpected missing pair for site at %s:%li. Skipping...", site->chr, site->pos+1);
            continue;
        }

        site->has_pair = 1;
        site->pair_position = pair.pos;
    }

    // Get competitors
    const int competitor_start = 2;
    const int competitor_end = sites.n-2;
    for (int i = 0; i < sites.n; i++) {
        SpliceSite *site = &sites.a[i];

        SpliceSite competitor;
        int competitor_missing;
        if (site->type == DONOR) {
            if (site->strand == STRAND_FWD) {
                if (i < competitor_start) continue;
                competitor_missing = find_next_site_with_type(sites, site->tid, i, competitor_start, -1, DONOR, &competitor);
            } else { // site.strand == STRAND_REV
                if (i >= competitor_end) continue;
                competitor_missing = find_next_site_with_type(sites, site->tid, i, competitor_end, 1, DONOR, &competitor);
            }
        } else {
            if (site->strand == STRAND_FWD) {
                if (i >= competitor_end) continue;
                competitor_missing = find_next_site_with_type(sites, site->tid, i, competitor_end, 1, ACCEPTOR, &competitor);
            } else { // site.strand == STRAND_REV
                if (i < competitor_start) continue;
                competitor_missing = find_next_site_with_type(sites, site->tid, i, competitor_start, -1, ACCEPTOR, &competitor);
            }
        }

        if (competitor_missing) {
            continue; // In contrast to missing pairs, missing competitors will occur in every exon around the start and end, so no warning
        }

        site->has_competitor = 1;
        site->competitor_position = competitor.pos;
    }

    return sites;
}

gff_t *read_gff(const char *gff_path) {
    log_debug("Loading GFF from: %s", gff_path);
    gff_t *gff = gff_init(gff_path);
    gff_parse(gff);
    return gff;
}


int predict_ref_position(const Model *models, const faidx_t *fai, const char *chr, const uint64_t position, const Range gene_bounds, const int strand, const int window_padding, int *out_num, float **out) {
    int window_len;
    const uint64_t window_beg = position - (SPLICEAI_CONTEXT_PADDING + window_padding);
    const uint64_t window_end = position + (SPLICEAI_CONTEXT_PADDING + window_padding);
    char *window_seq = faidx_fetch_seq(fai, chr, window_beg, window_end, &window_len);
    window_seq[window_len] = '\0';

    // Replace sequence beyond gene boundaries with N's
    const Range padding_size = {
        ((int64_t) gene_bounds.start - (int64_t) window_beg) > 0 ? (gene_bounds.start - window_beg) : 0,
        ((int64_t) window_end - (int64_t) gene_bounds.end) > 0 ? (window_end - gene_bounds.end) : 0
    };
    const char *padded_seq = pad_sequence(window_seq, padding_size, window_len);

    // Encode sequence
    float *encoded_seq = NULL;
    const int enc_len = one_hot_encode(padded_seq, window_len, (float **) &encoded_seq);
    if (strand == STRAND_REV) reverse_encoding(encoded_seq, enc_len);

    // Predict using encoding
    int num_predictions = 0;
    float *predictions = NULL;
    predict((Model *) models, enc_len, 1, &encoded_seq, &num_predictions, &predictions);
    free(encoded_seq);

    if (num_predictions != NUM_SCORES * (2*window_padding+1)) {
        log_error("Incorrect number of predictions made on site %s:%li, %i predictions instead of 3. This is a bug. Exiting program.", chr, position, num_predictions);
        return 1;
    }

    *out_num = num_predictions;
    *out = predictions;

    free(window_seq);

    return 0;
}

int predict_alt_position(const Model *models, const faidx_t *fai, const char *chr, const uint64_t position, const uint64_t variant_position, const char alt, const Range gene_bounds, const int strand, const int window_padding, int *out_num, float **out) {
    int window_len;
    const uint64_t window_beg = position - (SPLICEAI_CONTEXT_PADDING + window_padding);
    const uint64_t window_end = position + (SPLICEAI_CONTEXT_PADDING + window_padding);
    char *window_seq = faidx_fetch_seq(fai, chr, window_beg, window_end, &window_len);
    window_seq[window_len] = '\0';

    if ((window_beg > variant_position) | (window_end <= variant_position)) {
        return 1;
    }

    window_seq[variant_position - window_beg] = alt;

    // Replace sequence beyond gene boundaries with N's
    const Range padding_size = {
        ((int64_t) gene_bounds.start - (int64_t) window_beg) > 0 ? (gene_bounds.start - window_beg) : 0,
        ((int64_t) window_end - (int64_t) gene_bounds.end) > 0 ? (window_end - gene_bounds.end) : 0
    };
    const char *padded_seq = pad_sequence(window_seq, padding_size, window_len);

    // Encode sequence
    float *encoded_seq = NULL;
    const int enc_len = one_hot_encode(padded_seq, window_len, (float **) &encoded_seq);
    if (strand == STRAND_REV) reverse_encoding(encoded_seq, enc_len);

    // Predict using encoding
    int num_predictions = 0;
    float *predictions = NULL;
    predict((Model *) models, enc_len, 1, &encoded_seq, &num_predictions, &predictions);
    free(encoded_seq);

    if (num_predictions != NUM_SCORES * (2*window_padding+1)) {
        log_error("Incorrect number of predictions made on site %s:%li, %i predictions instead of 3. This is a bug. Exiting program.", chr, position, num_predictions);
        return 1;
    }

    *out_num = num_predictions;
    *out = predictions;

    free(window_seq);

    return 0;
}

void get_diff_position(const size_t num_preds, const float *refs, const float *alts, const SpliceSiteType type, int *out_index) {
    int max_gain_index = 0;
    float max_gain = 0.0;
    for (int i = type; i < num_preds; i += NUM_SCORES) {
        const float gain = alts[i] - refs[i];

        if (gain <= max_gain) continue;

        max_gain = gain;

        max_gain_index = ((i - type) - (SPLICEAI_WINDOW_PADDING * NUM_SCORES)) / NUM_SCORES; // transform to an index of -50:50
    }

    *out_index = max_gain_index;
}

CompactOutputRow init_compact_row(
    const char *chr, const uint64_t pos, const char *gene, const int tid, const int strand, const SpliceSiteType type,
    const uint64_t var_pos, const float var_ref, const float var_alt, const int var_dp,
    const float canonical_ref, const float canonical_alt,
    const int has_pair, const float pair_ref, const float pair_alt,
    const int has_competitor, const float competitor_ref, const float competitor_alt
) {
    return (CompactOutputRow) { chr, pos, gene, tid, strand, type, var_pos, var_ref, var_alt, var_dp, canonical_ref, canonical_alt, has_pair, pair_ref, pair_alt, has_competitor, competitor_ref, competitor_alt };
}

void build_compact_output_line(CompactOutputRow row, kstring_t *s) {
    kputs(row.chr, s);
    kputc('\t', s);
    kputl(row.pos+1, s);
    kputc('\t', s);
    kputs(row.gene_name, s);
    kputc('\t', s);
    kputw(row.tid, s);
    kputc('\t', s);
    kputs(row.strand == STRAND_FWD ? "fwd" : "rev", s);
    kputc('\t', s);
    kputs(row.type == ACCEPTOR ? "acceptor" : "donor", s);
    kputc('\t', s);
    kputw(row.var_pos, s);
    kputc('\t', s);
    kputd(row.var_ref, s);
    kputc('\t', s);
    kputd(row.var_alt, s);
    kputc('\t', s);
    kputw(row.var_dp, s);
    kputc('\t', s);
    kputd(row.canonical_ref, s);
    kputc('\t', s);
    kputd(row.canonical_alt, s);
    kputc('\t', s);
    if (row.has_pair) kputd(row.pair_ref, s);
    kputc('\t', s);
    if (row.has_pair) kputd(row.pair_alt, s);
    kputc('\t', s);
    if (row.has_competitor) kputd(row.competitor_ref, s);
    kputc('\t', s);
    if (row.has_competitor) kputd(row.competitor_alt, s);
}

ProbabilityOutputRow init_probability_row(const char *chr, const uint64_t pos, const char *gene, const int tid, const int strand, const SpliceSiteType type, const ProbabilityType prob_type, float *probabilities) {
    return (ProbabilityOutputRow) { chr, pos, gene, tid, strand, type, prob_type, probabilities };

}

void build_probabilities_output_line(ProbabilityOutputRow row, kstring_t *s) {
    kputs(row.chr, s);
    kputc('\t', s);
    kputl(row.pos+1, s);
    kputc('\t', s);
    kputs(row.gene_name, s);
    kputc('\t', s);
    kputw(row.tid, s);
    kputc('\t', s);
    kputs(row.strand == STRAND_FWD ? "fwd" : "rev", s);
    kputc('\t', s);
    kputs(row.type == ACCEPTOR ? "acceptor" : "donor", s);
    kputc('\t', s);
    kputs(row.prob_type == REF ? "ref" : "alt", s);

    for (int i = 0; i <= 101; i++) {
        kputc('\t', s);
        kputd(row.probabilities[i * NUM_SCORES + row.type], s);
    }
}


int main(int argc, char *argv[]) {
    setenv("TF_CPP_MIN_LOG_LEVEL", "1", TENSORFLOW_LOG_LEVEL_SILENT);

    if (argc != 5) {
        fprintf(stdout, "Run as: ./canonical_splice_analyzer <spliceai_model_dir> <human_fa> <gff> <out.tsv>");
    }

    const char *models_dir_path = argv[1];
    log_debug("Loading SpliceAI models from: %s", models_dir_path);
    Model *models = load_models(argv[1]);

    const char *fasta_path = argv[2];
    log_debug("Loading FASTA from: %s", fasta_path);
    faidx_t *fai = fai_load(fasta_path);

    gff_t *gff = read_gff(argv[3]);
    SpliceSites sites = get_splice_sites_from_gff(gff);
    log_debug("Found %i splice sites", sites.n);

    FILE *out = fopen(argv[4], "w");
    fprintf(out, COMPACT_HEADER_STRING);

    FILE *out_scores = fopen(argv[5], "w");
    kstring_t header = {0};
    kputs(PREFIX_HEADER_STRING, &header);
    for (int i = 0; i <= 101; i++) {
        kputc('\t', &header);
        kputs("idx_", &header);
        kputw(i, &header);
    }
    fprintf(out_scores, "%s\n", header.s);

    for (int i = 0; i < sites.n; i++) {
        SpliceSite site = sites.a[i];
        const Range gene_bounds = (Range) {site.gene_beg, site.gene_end};

        // Ref predictions
        int num_ref_canonical_predictions;
        float *ref_canonical_predictions;
        predict_ref_position(models, fai, site.chr, site.pos, gene_bounds, site.strand, SPLICEAI_WINDOW_PADDING, &num_ref_canonical_predictions, &ref_canonical_predictions);

        int num_ref_pair_predictions;
        float *ref_pair_predictions;
        if (site.has_pair) {
            predict_ref_position(models, fai, site.chr, site.pair_position, gene_bounds, site.strand, 0, &num_ref_pair_predictions, &ref_pair_predictions);
        }

        int num_ref_competitor_predictions;
        float *ref_competitor_predictions;
        if (site.has_competitor) {
            predict_ref_position(models, fai, site.chr, site.competitor_position, gene_bounds, site.strand, 0, &num_ref_competitor_predictions, &ref_competitor_predictions);
        }

        // Get positions from splice site to mutate
        Range splice_site_range;
        if (site.type == DONOR) {
            if (site.strand == STRAND_FWD) {
                splice_site_range = (Range) { site.pos + 1, site.pos + 3 };
            } else { // site.strand == STRAND_REV
                splice_site_range = (Range) { site.pos - 1, site.pos + 1 };
            }
        } else { // site.type == ACCEPTOR
            if (site.strand == STRAND_FWD) {
                splice_site_range = (Range) { site.pos - 1, site.pos + 1};
            } else { // site.strand == STRAND_REV
                splice_site_range = (Range) { site.pos + 1, site.pos + 3 };
            }
        }

        // Retrieve reference bases for splice site
        int splice_site_len;
        char *site_seq = faidx_fetch_seq(fai, site.chr, splice_site_range.start, splice_site_range.end, &splice_site_len);
        site_seq[splice_site_len] = '\0';

        for (int p = 0; p < 2; p++) {
            const char ref = site_seq[p];
            for (int b = 0; b < 4; b++) {
                const char alt = BASES[b];
                if (ref == alt) continue;

                int num_alt_canonical_predictions;
                float *alt_canonical_predictions;
                predict_alt_position(models, fai, site.chr, site.pos, splice_site_range.start + p, alt, gene_bounds, site.strand, SPLICEAI_WINDOW_PADDING, &num_alt_canonical_predictions, &alt_canonical_predictions);

                int gain_position = 0;
                get_diff_position((2*SPLICEAI_WINDOW_PADDING+1) * NUM_SCORES, ref_canonical_predictions, alt_canonical_predictions, site.type, &gain_position);

                int gain_index = gain_position * NUM_SCORES + (NUM_SCORES * SPLICEAI_WINDOW_PADDING) + site.type;
                float ref_canonical_gain = ref_canonical_predictions[gain_index];
                float alt_canonical_gain = alt_canonical_predictions[gain_index];

                int num_alt_pair_predictions;
                float *alt_pair_predictions;
                if (site.has_pair) {
                    int variant_out_of_range = predict_alt_position(models, fai, site.chr, site.pair_position, splice_site_range.start + p, alt, gene_bounds, site.strand, 0, &num_alt_pair_predictions, &alt_pair_predictions);
                    if (variant_out_of_range) {
                        alt_pair_predictions = ref_pair_predictions;
                    }
                }

                int num_alt_competitor_predictions;
                float *alt_competitor_predictions;
                if (site.has_competitor) {
                    int variant_out_of_range = predict_alt_position(models, fai, site.chr, site.competitor_position, splice_site_range.start + p, alt, gene_bounds, site.strand, 0, &num_alt_competitor_predictions, &alt_competitor_predictions);
                    if (variant_out_of_range) {
                        alt_competitor_predictions = ref_competitor_predictions;
                    }
                }

                {
                    kstring_t s = {0};
                    float ref_pair = 0.0, alt_pair = 0.0, ref_competitor = 0.0, alt_competitor = 0.0;
                    if (site.has_pair) {
                        ref_pair = ref_pair_predictions[site.type];
                        alt_pair = alt_pair_predictions[site.type];
                    }
                    if (site.has_competitor) {
                        ref_competitor = ref_competitor_predictions[site.type];
                        alt_competitor = alt_competitor_predictions[site.type];
                    }
                    CompactOutputRow row = init_compact_row(
                        site.chr, site.pos, site.gene_name, site.tid, site.strand, site.type, splice_site_range.start + p, 
                        ref_canonical_gain, alt_canonical_gain, gain_position, 
                        ref_canonical_predictions[150 + site.type], alt_canonical_predictions[150 + site.type], 
                        site.has_pair, ref_pair, alt_pair,
                        site.has_competitor, ref_competitor, alt_competitor
                    );
                    build_compact_output_line(row, &s);
                    fprintf(out, "%s\n", s.s);
                }

                {
                    kstring_t s = {0};
                    ProbabilityOutputRow row = init_probability_row(site.chr, site.pos, site.gene_name, site.tid, site.strand, site.type, REF, ref_canonical_predictions);
                    build_probabilities_output_line(row, &s);
                    fprintf(out_scores, "%s\n", s.s);
                }
                {
                    kstring_t s = {0};
                    ProbabilityOutputRow row = init_probability_row(site.chr, site.pos, site.gene_name, site.tid, site.strand, site.type, ALT, alt_canonical_predictions);
                    build_probabilities_output_line(row, &s);
                    fprintf(out_scores, "%s\n", s.s);
                }
            }
        }

        free(site_seq);
    }

    fclose(out);
    fai_destroy(fai);
    splice_sites_destroy(&sites);
    gff_destroy(gff);

    return 0;
}

