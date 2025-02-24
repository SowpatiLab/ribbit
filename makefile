CXX      = g++	# GNU c++ compiler
CXXFLAGS = -O3	# optimisation level flag

BOOST_LIB  = -lboost_system		# include the boost library
SRC_SSW    = src/ssw.c src/ssw_cpp.cpp	# source for striped-smithwaterman alignment
# list of ribbit source files
SRC_RIBBIT = src/global_variables.cpp src/concatenate_output.cpp src/output_utils.cpp src/process_cigar.cpp src/parse_seed.cpp src/parse_smallmotif_seed.cpp src/merge_types.cpp src/parse_anchored_shiftxor.cpp src/parse_substitute_shiftxor.cpp src/parse_perfect_shiftxor.cpp src/bitseq_utils.cpp src/fasta_utils.cpp src/ribbit.cpp
BOOST_PROGRAM_OPTIONS_LIB = -lboost_program_options		# including the program options library from boost

# if there is a change in any of the ribbit source file make builds the executable
ribbit: $(SRC_RIBBIT)
	$(CXX) $(CXXFLAGS) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) $(BOOST_PROGRAM_OPTIONS_LIB) -o ribbit