g++ -O3 -std=c++1z -I/opt/homebrew/Cellar/boost/1.86.0_2/include/ -L/opt/homebrew/Cellar/boost/1.86.0_2/lib ssw.c ssw_cpp.cpp global_variables.cpp concatenate_output.cpp process_cigar.cpp parse_seed.cpp parse_smallmotif_seed.cpp \
                                  merge_types.cpp parse_anchored_shiftxor.cpp parse_substitute_shiftxor.cpp parse_perfect_shiftxor.cpp \
                                  bitseq_utils.cpp fasta_utils.cpp ribbit.cpp -o ribbit -lboost_program_options
