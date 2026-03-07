/// \file CaloXRunAction.cc
/// \brief Implementation of the CaloXRunAction class

#include "CaloXRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloXRunAction::CaloXRunAction(CaloXTree *histo)
    : G4UserRunAction(),
      hh(histo)
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloXRunAction::~CaloXRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloXRunAction::BeginOfRunAction(const G4Run *run)
{
  std::cout << "### Run " << run->GetRunID() << " start." << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloXRunAction::EndOfRunAction(const G4Run * /*aRun*/)
{
  // print histogram statistics
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
