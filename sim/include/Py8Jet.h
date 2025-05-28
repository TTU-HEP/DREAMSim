//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 20 16:45:59 2020 by ROOT version 6.20/04
// from TTree py8tree/Pythia8 Tree
// found on file: py8_qq_jets.root
//////////////////////////////////////////////////////////

#ifndef Py8Jet_h
#define Py8Jet_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

class Py8Jet
{
public:
   TTree *fChain;  //! pointer to the analyzed TTree or TChain
   Int_t fCurrent; //! current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   std::vector<int> *particle_ID;
   std::vector<double> *particle_E;
   std::vector<double> *particle_px;
   std::vector<double> *particle_py;
   std::vector<double> *particle_pz;

   // List of branches
   TBranch *b_particle_ID; //!
   TBranch *b_particle_px; //!
   TBranch *b_particle_py; //!
   TBranch *b_particle_pz; //!
   TBranch *b_particle_e;  //!

   Py8Jet(TTree *tree = 0);
   virtual ~Py8Jet();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop();
   virtual Bool_t Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef Py8Jet_cxx
Py8Jet::Py8Jet(TTree *tree) : fChain(0)
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   if (tree == 0)
   {
      TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("py8_qq_jets.root");
      if (!f || !f->IsOpen())
      {
         f = new TFile("py8_qq_jets.root");
      }
      f->GetObject("py8tree", tree);
   }
   Init(tree);
}

Py8Jet::~Py8Jet()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t Py8Jet::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t Py8Jet::LoadTree(Long64_t entry)
{
   // Set the environment to read one entry
   if (!fChain)
      return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0)
      return centry;
   if (fChain->GetTreeNumber() != fCurrent)
   {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Py8Jet::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   particle_ID = 0;
   particle_E = 0;
   particle_px = 0;
   particle_py = 0;
   particle_pz = 0;
   // Set branch addresses and branch pointers
   if (!tree)
      return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("particle_ID", &particle_ID, &b_particle_ID);
   fChain->SetBranchAddress("particle_px", &particle_px, &b_particle_px);
   fChain->SetBranchAddress("particle_py", &particle_py, &b_particle_py);
   fChain->SetBranchAddress("particle_pz", &particle_pz, &b_particle_pz);
   fChain->SetBranchAddress("particle_E", &particle_E, &b_particle_e);
   Notify();
}

Bool_t Py8Jet::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Py8Jet::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t Py8Jet::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}
#endif // #ifdef Py8Jet_cxx
