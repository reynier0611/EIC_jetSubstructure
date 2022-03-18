#include<iostream>
using std::cout;
using std::cerr;
using std::endl;

#include<string>
using std::string;

#include<vector>
using std::vector;

//#include "fastjet/ClusterSequence.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "treeUtilities.h"

//#include "testCanyonlandsDIS.h"
//#include "testDeathvalleyDIS.h"


int main(int argc, char* argv[])
{
  if(argc != 2)
    {
      cerr << "Wrong number of arguments" << endl;
      exit(EXIT_FAILURE);
    }

  //const int run = atoi(argv[1]);
  const char* output = argv[1];

  //const char* input="blah.root";
  //const char* input="/gpfs02/eic/bpage/home/ATHENA/fullSimuFiles/canyonlands-2.1/pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_5_096.0002.root";
  //const char* output="test.hist.root";

  //TString inputFiles=Form("/gpfs02/eic/bpage/home/ATHENA/fullSimuFiles/canyonlands-2.1/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d_*.root",run);
  //TString inputFiles=Form("/gpfs02/eic/bpage/home/ATHENA/fullSimuFiles/deathvalley-1.0/pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d_*.root",run);
  //TString inputFiles="/gpfs02/eic/sjoosten/pythia8NCDIS_18x275_minQ2_1_hiDiv_1-birksfix.root";

  //cout << "Input Files = " << inputFiles << endl;
  cout << "Output File = " << output << endl;

  TChain *myChain = new TChain("eventT");
  //myChain->Add(Form("/gpfs02/eic/bpage/home/ATHENA/fullSimuFiles/canyonlands-2.1/pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_%d_*.root",));
  myChain->Add("/direct/eic+u/bpage/ATHENA/jetWork/generatorSubstructure/genSimuTrees/test.hist.root");

  //cout << "Added File to Tree" << endl;

  // Instatiate Class
  //testCanyonlandsDIS t;
  //testCanyonlandsDIS t(myChain);
  //testDeathvalleyDIS t(myChain);
  disEventTree etr;
  etr.initRead(myChain);

  // Open Output Root File
  TFile *ofile = TFile::Open(output,"recreate");

  // Histograms
  TH1D *numParts = new TH1D("numParts","",200,0.,200.);
  TH1D *testPx = new TH1D("testPx","",500,-50.,50.);


  // Entries
  Long64_t nentries = etr.fChain->GetEntries();
  cout << nentries << endl;

  Long64_t nb = 0;
  for(Long64_t jentry=0; jentry<nentries; jentry++)
    {
      Long64_t ientry = etr.LoadTree(jentry);
      if(ientry < 0) break;

      nb = etr.GetEntry(jentry);

      numParts->Fill(etr.part_);

      for(int i=0; i<etr.part_; i++)
	{
	  testPx->Fill(etr.part_px[i]);
	}

    }

  // Write and Close Root File
  ofile->Write();
  ofile->Close();

  return 0;
}
