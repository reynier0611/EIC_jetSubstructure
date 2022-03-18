
ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = c++
endif

CXXFLAGS     = -std=c++11
CXXFLAGS     += `root-config --cflags`
#CXXFLAGS     += `fastjet-config --cxxflags`
#CXXFLAGS     += -I$(PYTHIA8)/include

LDFLAGS      = 
LDFLAGS     += `root-config --libs`
#LDFLAGS     += `fastjet-config --libs --plugins`
#LDFLAGS     += -L$(PYTHIA8)/lib -lpythia8

IDIR = ../
DEPS = treeUtilities.h

all : run.exe

#run.exe: analyzeJetSubstructure.o 
#	$(CXX) analyzeJetSubstructure.o -o run.exe $(LDFLAGS) 

#analyzeJetSubstructure.o : analyzeJetSubstructure.cxx
#	$(CXX) $(CXXFLAGS) -c analyzeJetSubstructure.cxx -o analyzeJetSubstructure.o -I.

run.exe: analyzeJetAxes.o
	$(CXX) analyzeJetAxes.o -o run.exe $(LDFLAGS)

analyzeJetAxes.o: analyzeJetAxes.cxx $(IDIR)/$(DEPS)
	$(CXX) $(CXXFLAGS) -c analyzeJetAxes.cxx -I. -I..

clean :
	rm -vf *.o *.exe *~
