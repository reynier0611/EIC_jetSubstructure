// Histograms

#ifndef HISTOGRAMUTILITIES
#define HISTOGRAMUTILITIES

#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/Recluster.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"

using namespace std;
using namespace fastjet;

// Species Names
const char *speciesNames[]={"charged","gamma","neutral"};


// Particle Plots
struct partCollection {

  // Particle Distributions
  TH2D *partPtVsEta[3];
  TH2D *partPVsEta[3];

  void init() {

    for(int i=0; i<3; i++)
      {
	partPtVsEta[i] = new TH2D(Form("partPtVsEta_%s",speciesNames[i]),Form("Particle Pt Vs Eta: %s;Eta;Pt",speciesNames[i]),100,-5.,5.,200,0.,50.);
	partPVsEta[i] = new TH2D(Form("partPVsEta_%s",speciesNames[i]),Form("Particle P Vs Eta: %s;Eta;Momentum",speciesNames[i]),100,-5.,5.,200,0.,50.);
      }
  }

  void fillParticle(TLorentzVector p, Int_t pdg) {

    if(pdg == 211 || pdg == 321 || pdg == 2212)
      {
	partPtVsEta[0]->Fill(p.PseudoRapidity(),p.Perp());
	partPVsEta[0]->Fill(p.PseudoRapidity(),p.P());
      }

    if(pdg == 22)
      {
	partPtVsEta[1]->Fill(p.PseudoRapidity(),p.Perp());
	partPVsEta[1]->Fill(p.PseudoRapidity(),p.P());
      }

    if(pdg == 2112 || pdg == 130)
      {
	partPtVsEta[2]->Fill(p.PseudoRapidity(),p.Perp());
	partPVsEta[2]->Fill(p.PseudoRapidity(),p.P());
      }
  }

};


// Standard Jet Plots
struct jetCollection {

  // Jet Quantities
  TH1D *jetPtHist; 
  TH1D *jetEtaHist;
  TH2D *jetPtVsEtaHist;
  TH2D *jetEVsEtaHist;

  void init(const char *name, const char* title)
  {
    jetPtHist = new TH1D(Form("jetPt_%s",name),Form("Jet Pt: %s",title),200,0.,100.);
    jetEtaHist = new TH1D(Form("jetEta_%s",name),Form("Jet Eta: %s",title),100,-5.,5.);
    jetPtVsEtaHist = new TH2D(Form("jetPtVsEta_%s",name),Form("Jet Pt Vs Eta: %s",title),100,-5.,5.,200,0.,100.);
    jetEVsEtaHist = new TH2D(Form("jetEVsEta_%s",name),Form("Jet Energy Vs Eta: %s",title),100,-5.,5.,300,0.,300.);
  }

  void fillJetCollection(const PseudoJet jet)
  {
    jetPtHist->Fill(jet.pt());
    jetEtaHist->Fill(jet.eta());
    jetPtVsEtaHist->Fill(jet.eta(),jet.pt());
    jetEVsEtaHist->Fill(jet.eta(),jet.e());
  }
};


// Softdrop Jet Plots
struct jetSoftCollection {

  // Jet Quantities
  TH2D *jetGroomVsOrigPtHist;
  TH1D *jetDeltaRHist; 
  TH1D *jetZHist;
  TH2D *jetDeltaRVsZHist;
  TH2D *testParametersHist;

  void init(const char *name, const char* title)
  {
    jetGroomVsOrigPtHist = new TH2D(Form("jetGroomVsOrigPt_%s",name),Form("Groomed Vs Original Jet Pt: %s",title),200,0.,100.,200,0.,100.);
    jetDeltaRHist = new TH1D(Form("jetDeltaR_%s",name),Form("Delta R Between Subjets: %s",title),1000,0.,2.);
    jetZHist = new TH1D(Form("jetZ_%s",name),Form("Symmetry Measure of Subjets: %s",title),1000,0.,2.);
    jetDeltaRVsZHist = new TH2D(Form("jetDeltaRVsZ_%s",name),Form("Delta R Vs Symmetry Measure of Subjets: %s",title),1000,0.,2.,1000,0.,2.);
    testParametersHist = new TH2D(Form("testParameters_%s",name),Form("Z and Beta Test: %s",title),10,0.,0.5,10,0.,5.);
  }

  void fillSoftJetCollection(const PseudoJet jet, const PseudoJet sdJet, double z, double beta)
  {
    jetGroomVsOrigPtHist->Fill(jet.pt(),sdJet.pt());
    jetDeltaRHist->Fill(sdJet.structure_of<contrib::SoftDrop>().delta_R());
    jetZHist->Fill(sdJet.structure_of<contrib::SoftDrop>().symmetry());
    jetDeltaRVsZHist->Fill(sdJet.structure_of<contrib::SoftDrop>().symmetry(),sdJet.structure_of<contrib::SoftDrop>().delta_R());
    testParametersHist->Fill(z,beta);
  }
};

#endif
