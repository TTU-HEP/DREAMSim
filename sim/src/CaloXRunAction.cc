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
  fTimer.Start(); // start timing the run
}


void CaloXRunAction::EndOfRunAction(const G4Run *aRun)
{
  // print histogram statistics
  fTimer.Stop(); // stop timing the run

  G4int nEvents = aRun->GetNumberOfEvent();
  G4double realT = fTimer.GetRealElapsed();   // wall-clock time [s]
  G4double cpuT = fTimer.GetUserElapsed() + fTimer.GetSystemElapsed(); // CPU time [s]

  std::cout << "### Run " << aRun->GetRunID() << " finished:  "
            << nEvents << " events processed." << std::endl;
  std::cout << "### Timing:  real (wall) = " << realT << " s"
            << ",  CPU (user+sys) = " << cpuT << " s"
            << "  (user = " << fTimer.GetUserElapsed()
            << " s, sys = " << fTimer.GetSystemElapsed() << " s)" << std::endl;
  if (nEvents > 0)
  {
    std::cout << "### Per event:  real = " << realT / nEvents * 1000.0 << " ms"
              << ",  CPU = " << cpuT / nEvents * 1000.0 << " ms" << std::endl;
  }
}

