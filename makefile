CXX      = g++	# GNU c++ compiler
CXXFLAGS = -O3 -w # optimisation level flag; suppress warnings
SHARED_LIBS = -shared -fPIC $(shell python -m pybind11 --includes)	# shared library flags

BOOST_LIB     = -lboost_system
BOOST_AUXLIBS = -lboost_program_options -lboost_filesystem # including the program options library from boost
PTHREAD_LIB   = -lpthread
ZLIB	      = -lz # zlib library for gzip file handling

# library includes for striped-smithwaterman alignment
SRC_SSW    = src/ssw.c src/ssw_cpp.cpp

# list of ribbit source files
SRC_RIBBIT = src/global_variables.cpp src/concatenate_output.cpp src/binomial_thresholds.cpp src/cigar_utils.cpp src/output_utils.cpp \
             src/process_cigar.cpp src/parse_seed.cpp src/parse_smallmotif_seed.cpp src/largeseq_align.cpp src/merge_types.cpp \
			 src/parse_anchored_shiftxor.cpp src/parse_substitute_shiftxor.cpp src/parse_perfect_shiftxor.cpp src/seed_utils.cpp \
			 src/bitseq_utils.cpp src/fasta_utils.cpp

SRC_MAIN = src/ribbit.cpp

# Identify the operating system
OS := $(shell uname -s)
ifeq ($(OS),Darwin)
	CXXFLAGS      = -O3 -std=c++1z -w
	BOOST_VERSION = $(shell ls /opt/homebrew/Cellar/boost/ | tail -n 1)
	BOOST_LIB     = -L/opt/homebrew/Cellar/boost/$(BOOST_VERSION)/lib	# boost library path
	INCLUDE       = -I/opt/homebrew/Cellar/boost/$(BOOST_VERSION)/include/  # boost include path
	PROFILER_LIB  = -L$(shell brew --prefix gperftools)/lib  -I$(shell brew --prefix gperftools)/include -lprofiler
	BOOST_AUXLIBS = -lboost_program_options -lboost_filesystem		#i boost program options library path
endif

# if there is a change in any of the ribbit source file make builds the executable
ribbit: $(SRC_RIBBIT) $(SRC_MAIN)

ifeq ($(OS),Darwin)
	@echo "Operating System: macOS"
	@echo "Boost version identified: " ${BOOST_VERSION}
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) $(SRC_MAIN) $(BOOST_AUXLIBS) $(ZLIB) -o ribbit
else ifeq ($(OS),Linux)
	@echo "Operating System: Linux"
	$(CXX) $(CXXFLAGS) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) $(SRC_MAIN) $(BOOST_AUXLIBS) $(PTHREAD_LIB) $(ZLIB) -o ribbit
else
	@echo "Operating System: Unknown"
endif

pymodule: $(SRC_RIBBIT)
ifeq ($(OS),Darwin)
	@echo "Operating System: macOS"
	@echo "Boost version identified: " ${BOOST_VERSION}
	$(CXX) $(CXXFLAGS) $(SHARED_LIBS) $(INCLUDE) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) ./src/pyribbit.cpp $(BOOST_AUXLIBS) -o ribbit$(shell python3-config --extension-suffix) -undefined dynamic_lookup
else ifeq ($(OS),Linux)
	@echo "Operating System: Linux"
	$(CXX) $(CXXFLAGS) $(SHARED_LIBS) $(BOOST_LIB) $(SRC_SSW) $(SRC_RIBBIT) ./src/pyribbit.cpp $(BOOST_AUXLIBS) $(PTHREAD_LIB) -o ribbit$(shell python3-config --extension-suffix) -undefined dynamic_lookup
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

# clean removes the executable
clean:
	@rm -f ribbit
	@echo "Clean done!"

