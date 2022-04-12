// Generate Pythia Events and Output To a Flat Tree

#include "Pythia8/Pythia.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector2.h"

using namespace Pythia8;

#include "treeUtilities.h"


int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		cerr << "Wrong number of arguments" << endl;
		cerr << "program.exe steer out.hist.root" << endl;
		exit(EXIT_FAILURE);
	}

	// Set Parameters
	const char* rootOut = argv[2];
	//const char* rootOut="test.hist.root";

	// Check Parameters
	cout << "Root Output = " << rootOut << endl;

	// Open Root File
	TFile *ofile = TFile::Open(rootOut,"recreate");

	// Create Output Tree
	TTree *eventT = new TTree("eventT","Tree Holding Event and Particle Information");

	// Initilize Output Tree Structure
	disEventTree et;
	et.initWrite(eventT);


	// Histograms
	TH2D *phaseSpaceHist = new TH2D("phaseSpaceHist","",600,-5.,1.,400,0.,4.);
	TH1D *kinTestHist = new TH1D("kinTestHist","",1000,-1.,1.);
	TH1D *yHist = new TH1D("yHist","",600,-5.,1.);
	TH1D *w2Hist = new TH1D("w2Hist","",1000,0.,100.);
	TH1D *nuHist = new TH1D("nuHist","",1000,-5.,5.);
	TH2D *yVsNuHist = new TH2D("yVsNuHist","",1000,-5.,5.,600,-5.,1.);


	// Set Up Pythia
	Pythia8::Pythia p8;
	Pythia8::Event &event = p8.event;

	// Read in Steering File & Define Other Settings
	p8.readFile(argv[1]);
	//p8.readFile("dis_18x275.steer");
	p8.readString("Main:timesAllowErrors = 1000");

	// Initialize Pythia
	p8.init(); 


	// Run Event Generation
	int nevents = p8.mode("Main:numberOfEvents");
	for(int ev=0; ev<nevents; ev++)
	{
		if(!p8.next()) continue;

		// Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
		Pythia8::Vec4 pProton = event[1].p();
		Pythia8::Vec4 peIn    = event[4].p();
		Pythia8::Vec4 peOut   = event[6].p();
		Pythia8::Vec4 pPhoton = peIn - peOut;

		// Q2, W2, Bjorken x, y.
		double Q2    = - pPhoton.m2Calc();
		double W2    = (pProton + pPhoton).m2Calc();
		double x     = Q2 / (2. * pProton * pPhoton);
		double y     = (pProton * pPhoton) / (pProton * peIn);
		double nu    = pProton*pPhoton / event[1].m2();

		// Fill Test Histos
		phaseSpaceHist->Fill(std::log10(x),std::log10(Q2));
		kinTestHist->Fill(Q2 - 4.0*pProton.e()*peIn.e()*x*y);
		yHist->Fill(std::log10(y));
		w2Hist->Fill(TMath::Sqrt(W2));
		nuHist->Fill(std::log10(nu));
		yVsNuHist->Fill(std::log10(nu),std::log10(y));

		// Fill Tree
		et.code = p8.info.code();
		et.q2 = Q2;
		et.x = x;
		et.y = y;
		et.w2 = W2;
		et.nu = nu;
		et.part_ = p8.event.size();

		for(int i=0; i<p8.event.size(); i++)
		{
			et.part_index[i] = p8.event[i].index();
			et.part_id[i] = p8.event[i].id();
			et.part_status[i] = p8.event[i].status();
			et.part_mother1[i] = p8.event[i].mother1();
			et.part_mother2[i] = p8.event[i].mother2();
			et.part_px[i] = p8.event[i].px();
			et.part_py[i] = p8.event[i].py();
			et.part_pz[i] = p8.event[i].pz();
			et.part_e[i] = p8.event[i].e();
			et.part_mass[i] = p8.event[i].m();

			Int_t isFinal = -1;
			if(p8.event[i].isFinal()) isFinal = 1;
			else isFinal = 0;

			et.part_isFinal[i] = isFinal;
		}

		et.fill(eventT);
	}

	// List Statistics
	p8.stat();

	// Write and Close Root File
	ofile->Write();
	ofile->Close();

	return 0;
}
