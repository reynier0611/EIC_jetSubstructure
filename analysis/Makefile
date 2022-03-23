
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

all : run.exe

#run.exe: analyzeJetSubstructure.o 
#	$(CXX) analyzeJetSubstructure.o -o run.exe $(LDFLAGS) 

#analyzeJetSubstructure.o : analyzeJetSubstructure.cxx
#	$(CXX) $(CXXFLAGS) -c analyzeJetSubstructure.cxx -o analyzeJetSubstructure.o -I.

run.exe: analyzeSimuJetSubstructure.o
	$(CXX) analyzeSimuJetSubstructure.o -o run.exe $(LDFLAGS)

analyzeSimuJetSubstructure.o: analyzeSimuJetSubstructure.cxx $(IDIR)/$(DEPS) histogramUtilities.h
	$(CXX) $(CXXFLAGS) -c analyzeSimuJetSubstructure.cxx -I. -I..

clean :
	rm -vf *.o *.exe *~
