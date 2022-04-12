
ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = c++
endif

CXXFLAGS     = -std=c++11
CXXFLAGS     += `root-config --cflags`
CXXFLAGS     += `fastjet-config --cxxflags`
#CXXFLAGS     += -I$(PYTHIA8)/include

LDFLAGS      = 
LDFLAGS     += `root-config --libs`
LDFLAGS     += `fastjet-config --libs --plugins`
#LDFLAGS     += -L$(PYTHIA8)/lib -lpythia8
LDFLAGS     += -lRecursiveTools

IDIR = ../
DEPS = treeUtilities.h

all : code.exe

code.exe: analyzeSimuJetSubstructure.o
	$(CXX) analyzeSimuJetSubstructure.o -o code.exe $(LDFLAGS)

analyzeSimuJetSubstructure.o: analyzeSimuJetSubstructure.cxx $(IDIR)/$(DEPS) histogramUtilities.h
	$(CXX) $(CXXFLAGS) -c analyzeSimuJetSubstructure.cxx -I. -I..

clean :
	rm -vf *.o *.exe *~
