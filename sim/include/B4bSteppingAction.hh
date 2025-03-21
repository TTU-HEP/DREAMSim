//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4bSteppingAction.hh
/// \brief Definition of the B4bSteppingAction class

#ifndef B4bSteppingAction_h
#define B4bSteppingAction_h 1

#include "G4UserSteppingAction.hh"

#include "B4bEventAction.hh"
// class B4bEventAction;

class CaloID;
class CaloHit;
class CaloTree;

/// Stepping action class.
///
/// In UserSteppingAction() there are collected the energy deposit and track
/// lengths of charged particles in Absober and Gap layers and
/// updated in B4bRunData object.

class B4bSteppingAction : public G4UserSteppingAction
{
public:
  B4bSteppingAction(B4bEventAction *eventAction, CaloTree *histo);
  virtual ~B4bSteppingAction();

  virtual void UserSteppingAction(const G4Step *step);

private:
  B4bEventAction *fEventAction;
  CaloTree *hh;

  double getBirk(const G4Step *step);
  double getBirkHC(double dEStep, double step, double charge, double density);
  double getBirkL3(double dEStep, double step, double charge, double density);
  vector<double> UserCerenkov(const G4Step *step);

  //  SiPM PDE handling...
  void initPDE(int sipmType); // 1= J 6mm 6.0V,  2=J 6 mm 2.5V
  float getPDE(float lambda);

  int pdeLambdaMin, pdeLambdaMax;
  std::vector<float> sipmPDE;
  double findInvisible(const G4Step *step, bool verbose = false);
  void fillOPInfo(const G4Step *step, bool verbose = false);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
