PKGCONFIG	:= pkg-config
PYTHON		:= $(shell which python)
ROOTCONFIG	:= root-config

CXX		:= $(shell $(ROOTCONFIG) --cxx)
OPT		:= -g -Wall -std=c++11
CXXFLAGS	:= -c $(OPT)

# Boost flags
BOOSTLIBS    := -lboost_regex
BOOSTLIBS.CLI:= -lboost_program_options

# ROOT flags
ROOTCFLAGS	:= $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS	:= $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS	:= $(shell $(ROOTCONFIG) --glibs)

SRCS	:= $(wildcard *.cc)
TESTS	:= $(wildcard test_*.py)

# targets
$(SRCS:%.cc=%):%:	%.cc
	$(CXX) $(OPT) $(BOOSTLIBS.CLI) $(ROOTCFLAGS) $(ROOTLIBS) -o $@ $<

-include local.mk

$(TESTS:%.py=%):%:	%.py
	 $(PYTHON) -m unittest $(OPTS) $@
