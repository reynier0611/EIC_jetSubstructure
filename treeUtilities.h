// Define Flat Tree Structure and Helper Functions

#ifndef treeUtilities_h
#define treeUtilities_h

//#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>

class disEventTree {
 public:
  TTree *fChain; // Tree Being Read

  // Max Array Size
  static const Int_t MAXPARTS = 500;

  // Declare Leaf Types
  Int_t code;
  Float_t q2;
  Float_t x;
  Float_t y;
  Float_t w2;
  Float_t nu;
  Int_t part_;
  Int_t part_index[MAXPARTS];
  Int_t part_id[MAXPARTS];
  Int_t part_status[MAXPARTS];
  Int_t part_mother1[MAXPARTS];
  Int_t part_mother2[MAXPARTS];
  Int_t part_isFinal[MAXPARTS];
  Float_t part_px[MAXPARTS];
  Float_t part_py[MAXPARTS];
  Float_t part_pz[MAXPARTS];
  Float_t part_e[MAXPARTS];
  Float_t part_mass[MAXPARTS];
  

  // List Branches
  //TBranch *b_q2;
  //TBranch *b_x;
  //TBranch *b_part_;
  //TBranch *b_part_px;
  //TBranch *b_part_py;
  //TBranch *b_part_pz;
  //TBranch *b_part_e;
  //TBranch *b_part_isFinal;

  disEventTree();
  virtual ~disEventTree();
  void initWrite(TTree *myTtree);
  void fill(TTree *myTtree);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  void initRead(TTree *myTtree);
};

disEventTree::disEventTree()
{
  //initWrite(myTtree);
}


disEventTree::~disEventTree()
{

}


Int_t disEventTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}


Long64_t disEventTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  return centry;
}


void disEventTree::initWrite(TTree *myTtree)
{
  // Set Up Tree Struct For Writing

  myTtree->Branch("code",&code,"code/I");
  myTtree->Branch("q2",&q2,"q2/F");
  myTtree->Branch("x",&x,"x/F");
  myTtree->Branch("y",&y,"y/F");
  myTtree->Branch("w2",&w2,"w2/F");
  myTtree->Branch("nu",&nu,"nu/F");
  myTtree->Branch("part_",&part_,"part_/I");
  myTtree->Branch("part_index",part_index,"part_index[part_]/I");
  myTtree->Branch("part_id",part_id,"part_id[part_]/I");
  myTtree->Branch("part_status",part_status,"part_status[part_]/I");
  myTtree->Branch("part_mother1",part_mother1,"part_mother1[part_]/I");
  myTtree->Branch("part_mother2",part_mother2,"part_mother2[part_]/I");
  myTtree->Branch("part_isFinal",part_isFinal,"part_isFinal[part_]/I");
  myTtree->Branch("part_px",part_px,"part_px[part_]/F");
  myTtree->Branch("part_py",part_py,"part_py[part_]/F");
  myTtree->Branch("part_pz",part_pz,"part_pz[part_]/F");
  myTtree->Branch("part_e",part_e,"part_e[part_]/F");
  myTtree->Branch("part_mass",part_mass,"part_mass[part_]/F");
}


void disEventTree::initRead(TTree *myTtree)
{
  fChain = myTtree;

  fChain->SetBranchAddress("code",&code);
  fChain->SetBranchAddress("q2",&q2);
  fChain->SetBranchAddress("x",&x);
  fChain->SetBranchAddress("y",&y);
  fChain->SetBranchAddress("w2",&w2);
  fChain->SetBranchAddress("nu",&nu);
  fChain->SetBranchAddress("part_",&part_);
  fChain->SetBranchAddress("part_index",part_index);
  fChain->SetBranchAddress("part_id",part_id);
  fChain->SetBranchAddress("part_status",part_status);
  fChain->SetBranchAddress("part_mother1",part_mother1);
  fChain->SetBranchAddress("part_mother2",part_mother2);
  fChain->SetBranchAddress("part_isFinal",part_isFinal);
  fChain->SetBranchAddress("part_px",part_px);
  fChain->SetBranchAddress("part_py",part_py);
  fChain->SetBranchAddress("part_pz",part_pz);
  fChain->SetBranchAddress("part_e",part_e);
  fChain->SetBranchAddress("part_mass",part_mass);
}


void disEventTree::fill(TTree *myTtree)
{
  myTtree->Fill();
}

#endif
