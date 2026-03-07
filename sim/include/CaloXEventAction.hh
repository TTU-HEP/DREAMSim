/// \file CaloXEventAction.hh
/// \brief Definition of the CaloXEventAction class

#ifndef CaloXEventAction_h
#define CaloXEventAction_h 1

#include "G4UserEventAction.hh"
#include "CaloXPrimaryGeneratorAction.hh"
#include "CaloXDetectorConstruction.hh"
#include "G4SteppingManager.hh"

#include "globals.hh"

#include "G4Step.hh"

#include <stdlib.h> /* getenv */

using namespace std;

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy
/// deposit and track lengths of charged particles in Absober and Gap layers
/// stored in CaloXRunData object.

class CaloXTree;

class CaloXEventAction : public G4UserEventAction
{
public:
  CaloXEventAction(CaloXDetectorConstruction *det, CaloXPrimaryGeneratorAction *prim, CaloXTree *);
  virtual ~CaloXEventAction();

  virtual void BeginOfEventAction(const G4Event *event);
  virtual void EndOfEventAction(const G4Event *event);

private:
  // methods
  CaloXDetectorConstruction *fDetector;
  CaloXPrimaryGeneratorAction *primary;
  CaloXTree *hh;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
