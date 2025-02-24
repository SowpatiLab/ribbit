##
## Makefile for all executables
##

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXX      = g++
CXXFLAGS = -O3 -std=c++1z

SRC_SSW    = src/ssw.c src/ssw_cpp.cpp
SRC_RIBBIT = src/global_variables.cpp src/concatenate_output.cpp src/process_cigar.cpp src/parse_seed.cpp src/parse_smallmotif_seed.cpp src/merge_types.cpp src/parse_anchored_shiftxor.cpp src/parse_substitute_shiftxor.cpp src/parse_perfect_shiftxor.cpp src/bitseq_utils.cpp src/fasta_utils.cpp src/ribbit.cpp
INCLUDE    = -I/opt/homebrew/Cellar/boost/1.86.0_2/include/
BOOST_LIB  = -L/opt/homebrew/Cellar/boost/1.86.0_2/lib
BOOST_PROGRAM_OPTIONS_LIB = -lboost_program_options

ribbit: $(SRC_RIBBIT)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) $(BOOST_PROGRAM_OPTIONS_LIB) -o ribbit