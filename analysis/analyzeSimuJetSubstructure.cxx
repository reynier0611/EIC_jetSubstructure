#include<iostream>
using std::cout;
using std::cerr;
using std::endl;

#include<string>
using std::string;

#include<vector>
using std::vector;

#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "treeUtilities.h"
#include "histogramUtilities.h"

using namespace std;
using namespace fastjet;

// =======================================================================================================================================
// Main function
int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		cerr << "***************************************************" << endl;
		cerr << "Wrong number of arguments" << endl;
		cerr << "Run this code as:" << endl;
		cerr << "./analyzer /path/to/input.root /path/to/output.root" << endl;
		cerr << "***************************************************" << endl;
		exit(EXIT_FAILURE);
	}
	// ----------------------------------------------------------------------------------
	// Gathering some needed variables
	const char* inputFiles = argv[1];
	const char* output = argv[2];
	cout << "Input Files = " << inputFiles << endl;
	cout << "Output File = " << output << endl;
	
	// Jets
	double jetR[] = {1.0};
	const int size_jetR = sizeof(jetR)/sizeof(*jetR);
	double min_jet_pT = 1.0; // GeV
	// ----------------------------------------------------------------------------------

	TChain *myChain = new TChain("eventT");
	myChain->Add(inputFiles);

	// Instantiate Class
	disEventTree etr;
	etr.initRead(myChain);

	// Open Output Root File
	TFile *ofile = TFile::Open(output,"recreate");

	// Histograms
	TH1::SetDefaultSumw2(true);
	TH2::SetDefaultSumw2(true);

	// Event Kinematics
	TH2D *phaseSpaceHist = new TH2D("phaseSpace","Event Q2 Vs x",600,-5.,1.,400,0.,4.);
	TH1D *inelasticityHist = new TH1D("inelasticity","Event y",600,-5.,1.);
	TH1D *w2Hist = new TH1D("w2","Event W^2",2000,0.,200.);
	TH1D *nuHist = new TH1D("nu","Event nu",1000,-5.,5.);

	// Generated Particles
	partCollection pC;
	pC.init();

	// Jets
	char jetNames[2][20];
	sprintf(jetNames[0],"full");
	sprintf(jetNames[1],"track");

	jetCollection jC[2];
	for(int i=0; i<2; i++)
	{
		TString name = Form("%s",jetNames[i]);
		TString title = Form("%s",jetNames[i]);
		jC[i].init(name,title);
	}

	char jetZNames[3][20];
	sprintf(jetZNames[0],"z005");
	sprintf(jetZNames[1],"z01");
	sprintf(jetZNames[2],"z02");

	char jetBetaNames[3][20];
	sprintf(jetBetaNames[0],"beta0");
	sprintf(jetBetaNames[1],"beta1");
	sprintf(jetBetaNames[2],"beta2");

	jetSoftCollection jSC[2][9];
	for(int i=0; i<2; i++)
	{
		for(int j=0; j<3; j++)
		{
			for(int k=0; k<3; k++)
			{
				int flatIndex = 3*j+k;

				TString name = Form("%s_%s_%s",jetNames[i],jetZNames[j],jetBetaNames[k]);
				TString title = Form("%s %s %s",jetNames[i],jetZNames[j],jetBetaNames[k]);
				jSC[i][flatIndex].init(name,title);
			}
		}
	}


	TH1D *numParts = new TH1D("numParts","",200,0.,200.);
	TH1D *numPartsJet = new TH1D("numPartsJet","",200,0.,200.);


	// ----------------------------------------------------------------------------------
	// Loop Over Events
	Long64_t nentries = etr.fChain->GetEntries();

	Long64_t nb = 0;
	for(Long64_t jentry=0; jentry<nentries; jentry++)
	{
		// Load Tree and Get Entry
		Long64_t ientry = etr.LoadTree(jentry);
		if(ientry < 0) break;
		nb = etr.GetEntry(jentry);

		if(ientry%50000 == 0) cout << "Event " << ientry << endl;

		// Event Level Quantities
		phaseSpaceHist->Fill(etr.x,etr.q2);
		inelasticityHist->Fill(etr.y);
		w2Hist->Fill(etr.w2);
		nuHist->Fill(etr.nu);

		// Create FastJet Particle Containers
		vector<PseudoJet> particles;
		vector<PseudoJet> particlesTrack;

		numParts->Fill(etr.part_);
	
		// -------------------------------------------------------------
		// Loop Over Particles
		int numJetParts = 0;
		for(int i=0; i<etr.part_; i++)
		{
			// Fill Particle Level Plots
			Float_t px = etr.part_px[i];
			Float_t py = etr.part_py[i];
			Float_t pz = etr.part_pz[i];
			Float_t e = etr.part_e[i];
			Float_t m = etr.part_mass[i];
			Int_t final = etr.part_isFinal[i];
			Int_t pdg = etr.part_id[i];
			Int_t index = etr.part_index[i];

			// Particle Momentum
			TLorentzVector partP;
			partP.SetPxPyPzE(px,py,pz,e);

			// Populate FastJet
			if(final == 1 && TMath::Abs(partP.PseudoRapidity()) < 4.0)
			{
				if(pdg != 11)
				{
					fastjet::PseudoJet p(px,py,pz,e);
					p.set_user_index(index);
					particles.push_back(p);

					numJetParts++;

					if(TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 2212)
					{
						fastjet::PseudoJet pT(px,py,pz,e);
						pT.set_user_index(index);
						particlesTrack.push_back(pT);
					}

					pC.fillParticle(partP,pdg);
				}
			}
		}
		numPartsJet->Fill(numJetParts);

		// Find Jets
		// Loop Over Jet Radii
		for(int rad=0; rad<size_jetR; rad++)
		{
			// Define Jet
			JetDefinition jet_def_akt(antikt_algorithm,jetR[rad]);

			// Cluster Particles
			ClusterSequence cs_akt_lab(particles, jet_def_akt);
			ClusterSequence cs_akt_lab_track(particlesTrack, jet_def_akt);

			// Get Jets
			vector<vector<PseudoJet> > jetsVec; // Hold The Jet Vectors

			// Full jet
			vector<PseudoJet> jets_akt_lab = sorted_by_pt(cs_akt_lab.inclusive_jets(min_jet_pT));
			jetsVec.push_back(jets_akt_lab);

			// Charged jets
			vector<PseudoJet> jets_akt_lab_track = sorted_by_pt(cs_akt_lab_track.inclusive_jets(min_jet_pT));
			jetsVec.push_back(jets_akt_lab_track);

			// Loop Over Jet Types
			for(int jV=0; jV<jetsVec.size(); jV++)
			{
				// Loop over jets of a given type
				for(unsigned int jn=0; jn<jetsVec[jV].size(); jn++)
				{
					// Place Kin Cuts on Individual Jets
					if(jetsVec[jV][jn].pt() < 5.0) continue;
					if(TMath::Abs(jetsVec[jV][jn].eta()) > 3.0) continue;

					// Get Constituents for Recluster
					vector<PseudoJet> constituents = jetsVec[jV][jn].constituents();
					if(constituents.size() == 1) continue;

					// Fill Baseline Jet Histograms
					jC[jV].fillJetCollection(jetsVec[jV][jn]);

					// Recluster with Cambridge Aachen for Grooming
					double R25 = 2.5;
					JetDefinition jet_def_CA(fastjet::cambridge_algorithm,R25);
					ClusterSequence cs_CA(constituents, jet_def_CA);
					vector<PseudoJet> jets_CA = sorted_by_pt(cs_CA.inclusive_jets(min_jet_pT));

					// Setup Softdrop
					double zVals[3] = {0.05,0.1,0.2};
					double bVals[3] = {0.0,1.0,2.0};
					for(int zcut=0; zcut<3; zcut++)
					{
						for(int beta=0; beta<3; beta++)
						{
							contrib::SoftDrop sd( bVals[beta], zVals[zcut]);
							PseudoJet sd_jet = sd( jets_CA[0]);

							// Find Flat Index (Should be the same as in container initilization)
							int flatIndex = 3*zcut+beta;

							// Fill Soft-Drop Jet Histograms
							jSC[jV][flatIndex].fillSoftJetCollection(jetsVec[jV][jn],sd_jet,zVals[zcut],bVals[beta]);
						}
					}			  
				} // End Loop Through Jets
			} // End Loop Through Jet Types

		} // End Loop Over Jet Radii

	}

	// Write and Close Root File
	ofile->Write();
	ofile->Close();

	return 0;
}
// =======================================================================================================================================
