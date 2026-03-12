/// \file CaloXRunAction.cc
/// \brief Implementation of the CaloXRunAction class

#include "CaloXRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


CaloXRunAction::CaloXRunAction(CaloXTree *histo)
    : G4UserRunAction(),
      hh(histo)
{
  // Print progress every N events (0 = quiet). Increase for verbose output.
  G4RunManager::GetRunManager()->SetPrintProgress(100);
}


CaloXRunAction::~CaloXRunAction()
{
}


void CaloXRunAction::BeginOfRunAction(const G4Run *run)
{
  std::cout << "### Run " << run->GetRunID() << " start." << std::endl;
}


void CaloXRunAction::EndOfRunAction(const G4Run * /*aRun*/)
{
  // print histogram statistics
}

