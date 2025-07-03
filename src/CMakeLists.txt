set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -DLOG_USE_COLOR")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -DLOG_USE_COLOR")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

# Create executable
add_executable(splice_attenuation_analysis splice_attenuation_analysis.c)

find_package(Sanitizers)

# Link external libraries
target_link_libraries(splice_attenuation_analysis ${HTSLIB_LIBRARIES})
target_link_libraries(splice_attenuation_analysis ${CPLICEAI_LIBRARIES})

# Add logging
add_library(Logging STATIC logging/log.c)
target_compile_definitions(Logging PRIVATE -DLOG_USE_COLOR)
target_include_directories(splice_attenuation_analysis PRIVATE Logging)
target_link_libraries(splice_attenuation_analysis Logging)

# Add BCFTools
add_library(Bcftools STATIC bcftools/regidx.c bcftools/gff.c)
target_link_libraries(Bcftools PRIVATE ${HTSLIB_LIBRARIES})
target_include_directories(splice_attenuation_analysis PRIVATE Bcftools)
target_link_libraries(splice_attenuation_analysis Bcftools)

