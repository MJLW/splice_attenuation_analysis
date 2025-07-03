# canonical_splice_analyzer
Runs CpliceAI on all protein coding canonical splice sites
Requires CpliceAI and HTSLIB. Build using CMake:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..
cmake --build build
```


Run as:
```
./build/canonical_splice_analyzer <spliceai_model_dir> <human_fa> <gff> <out.tsv>
```
