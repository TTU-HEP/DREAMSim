/// \file CaloXPrimaryGeneratorAction.hh
/// \brief Definition of the CaloXPrimaryGeneratorAction class

#ifndef CaloXPrimaryGeneratorAction_h
#define CaloXPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "G4ParticleGun.hh"
#include "CaloXDetectorConstruction.hh"

#include <stdlib.h> /* getenv */

class TFile;
class Py8Jet;

class G4ParticleGun;
class G4Event;

class CaloXTree;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

class CaloXPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  CaloXPrimaryGeneratorAction(CaloXDetectorConstruction *det, CaloXTree *histo);
  virtual ~CaloXPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event *event);

  // set methods
  void SetRandomFlag(G4bool value);
  G4ParticleGun *GetParticleGun() { return fParticleGun; };

  void FillHEPparticles(std::vector<int> *mHepPID, std::vector<int> *mHepStatus,
                        std::vector<int> *mHepMother1, std::vector<int> *mHepMother2,
                        std::vector<int> *mHepDaughter1, std::vector<int> *mHepDaughter2,
                        std::vector<float> *mHepPx, std::vector<float> *mHepPy,
                        std::vector<float> *mHepPz, std::vector<float> *mHepE,
                        std::vector<float> *mHepMass);

private:
  G4ParticleGun *fParticleGun; // G4 particle gun
  CaloXDetectorConstruction *fDetector;
  CaloXTree *hh;

  G4ParticleTable *particleTable;

  double worldZHalfLength;
  void getPy8Event(G4Event *event);
  void printPy8Event();

  // parameters from env. variables
  void getParamFromEnvVars();
  int CaloXPythiaON;                       //  0=singl particle gun, 1=Pythia8 Root file.
  double CaloXPythiaXmin, CaloXPythiaXmax; //  vertex point smearing.
  double CaloXPythiaYmin, CaloXPythiaYmax;
  double CaloXPythiaZmin, CaloXPythiaZmax;
  int CaloXPythiaSkip;  //  number of events to skip
  int CaloXPythiaPrint; //  number of events to print
  std::string CaloXPythiaFile;

  TFile *finPy8;
  Py8Jet *py8evt;
  int py8eventCounter;
  int py8eventNumber;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
