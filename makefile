CXX      = g++	# GNU c++ compiler
CXXFLAGS = -O3 -Wall -g # optimisation level flag; suppress warnings
SHARED_LIBS = -Wall -shared -fPIC $(shell python -m pybind11 --includes)	# shared library flags

BOOST_LIB  = -lboost_system
BOOST_AUXLIBS = -lboost_program_options -lboost_filesystem # including the program options library from boost
PTHREAD_LIB = -lpthread

# library includes for striped-smithwaterman alignment
SRC_SSW    = test/ssw.c test/ssw_cpp.cpp

# list of ribbit source files
SRC_RIBBIT = test/global_variables.cpp test/concatenate_output.cpp test/binomial_thresholds.cpp test/cigar_utils.cpp test/output_utils.cpp \
             test/process_cigar.cpp test/parse_seed.cpp test/parse_smallmotif_seed.cpp test/merge_types.cpp \
			 test/parse_anchored_shiftxor.cpp test/parse_substitute_shiftxor.cpp test/parse_perfect_shiftxor.cpp test/seed_utils.cpp \
			 test/bitseq_utils.cpp test/fasta_utils.cpp

SRC_MAIN = test/ribbit.cpp
SRC_COMPLEX = test/complex_utils.cpp

# Identify the operating system
OS := $(shell uname -s)
ifeq ($(OS),Darwin)
	CXXFLAGS      = -O3 -std=c++1z -w -lz
	BOOST_VERSION = $(shell ls /opt/homebrew/Cellar/boost/ | tail -n 1)
	BOOST_LIB     = -L/opt/homebrew/Cellar/boost/$(BOOST_VERSION)/lib	# boost library path
	INCLUDE       = -I/opt/homebrew/Cellar/boost/$(BOOST_VERSION)/include/  # boost include path
	PROFILER_LIB  = -L$(shell brew --prefix gperftools)/lib  -I$(shell brew --prefix gperftools)/include -lprofiler
	BOOST_AUXLIBS = -lboost_program_options -lboost_filesystem		#i boost program options library path
endif

# if there is a change in any of the ribbit source file make builds the executable
ribbit: $(SRC_RIBBIT)

ifeq ($(OS),Darwin)
	@echo "Operating System: macOS"
	@echo "Boost version identified: " ${BOOST_VERSION}
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) $(SRC_MAIN) $(BOOST_AUXLIBS) -o ribbit
else ifeq ($(OS),Linux)
	@echo "Operating System: Linux"
	$(CXX) $(CXXFLAGS) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) $(SRC_MAIN) $(BOOST_AUXLIBS) $(PTHREAD_LIB) -o ribbit
else
	@echo "Operating System: Unknown"
endif

profile: $(SRC_RIBBIT)

ifeq ($(OS),Darwin)
	@echo "Operating System: macOS"
	@echo "Boost version identified: " ${BOOST_VERSION}
	$(CXX) -g -std=c++1z -w $(INCLUDE) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) $(SRC_MAIN) $(BOOST_AUXLIBS) -o ribbit $(PROFILER_LIB)
else ifeq ($(OS),Linux)
	@echo "Operating System: Linux"
	$(CXX) $(CXXFLAGS) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) $(SRC_MAIN) $(BOOST_AUXLIBS) $(PTHREAD_LIB) -o ribbit -lprofiler
else
	@echo "Operating System: Unknown"
endif

pymodule: $(SRC_RIBBIT)
ifeq ($(OS),Darwin)
	@echo "Operating System: macOS"
	@echo "Boost version identified: " ${BOOST_VERSION}
	$(CXX) $(CXXFLAGS) $(SHARED_LIBS) $(INCLUDE) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) ./test/pyribbit.cpp $(BOOST_AUXLIBS) -o ribbit$(shell python3-config --extension-suffix) -undefined dynamic_lookup
else ifeq ($(OS),Linux)
	@echo "Operating System: Linux"
	$(CXX) $(CXXFLAGS) $(SHARED_LIBS) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) ./test/pyribbit.cpp $(BOOST_AUXLIBS) $(PTHREAD_LIB) -o ribbit$(shell python3-config --extension-suffix) -undefined dynamic_lookup
else
	@echo "Operating System: Unknown"
endif


complex: $(SRC_COMPLEX)
ifeq ($(OS),Darwin)
	@echo "Operating System: macOS"
	@echo "Boost version identified: " ${BOOST_VERSION}
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(BOOST_LIB) $(SRC_SSW) $(SRC_COMPLEX) $(BOOST_AUXLIBS) -o ribbit_complex
else ifeq ($(OS),Linux)
	@echo "Operating System: Linux"
	$(CXX) $(CXXFLAGS) $(BOOST_LIB) $(SRC_SSW) $(SRC_COMPLEX) $(BOOST_AUXLIBS) $(PTHREAD_LIB) -o ribbit_complex
else
	@echo "Operating System: Unknown"
endif

# clean removes the executable
clean:
	@rm -f ribbit
	@echo "Clean done!"

