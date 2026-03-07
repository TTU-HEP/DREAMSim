/// \file CaloXEventAction.cc
/// \brief Implementation of the CaloXEventAction class

#include "CaloXEventAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include "G4Step.hh"

#include "CaloXPrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include <iomanip>

// -- for CaloX data --
#include "CaloXID.h"
#include "CaloXHit.h"
#include "CaloXTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloXEventAction::CaloXEventAction(CaloXDetectorConstruction *det, CaloXPrimaryGeneratorAction *prim, CaloXTree *histo)
    : G4UserEventAction(), fDetector(det), primary(prim), hh(histo)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloXEventAction::~CaloXEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloXEventAction::BeginOfEventAction(const G4Event *event /*event*/)
{
   // std::cout<<"CaloXEventAction::BeginOfEventAction-  starting..."<<std::endl;
   // G4Random::showEngineStatus();
   hh->BeginEvent();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloXEventAction::EndOfEventAction(const G4Event *event)
{
   // std::cout<<"CaloXEventAction::EndOfEventAction-  starting..."<<std::endl;
   hh->EndEvent();

} //  end of CaloXEventAction::EndOfEventAction

// -----------------------------------------------------------------------
