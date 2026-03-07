/// \file CaloXPrimaryGeneratorAction.cc
/// \brief Implementation of the CaloXPrimaryGeneratorAction class

#include "CaloXPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "Py8Jet.h"

#include "CaloXDetectorConstruction.hh"

// for root tree
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
// #include "TROOT.h"


#include "CaloXTree.h"

using namespace std;

CaloXPrimaryGeneratorAction::CaloXPrimaryGeneratorAction(CaloXDetectorConstruction *det, CaloXTree *histo)
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(nullptr), fDetector(det), hh(histo)
{
  cout << "CaloXPrimaryGeneratorAction constructer is called..." << endl;
  // Create the table containing all particle names
  particleTable = G4ParticleTable::GetParticleTable();

  getParamFromEnvVars(); // get paramters from theenvvariables.

  if (CaloXPythiaON == 1)
  {
    py8eventCounter = 0;
    py8eventNumber = CaloXPythiaSkip;
    string inFileName = CaloXPythiaFile;
    std::cout << "CaloXPrimaryGeneratorAction:  Using Pythia Event file: " << inFileName << std::endl;
    finPy8 = new TFile(inFileName.c_str());
    TTree *py8tree;
    finPy8->GetObject("Particles", py8tree);
    py8evt = new Py8Jet(py8tree);
    py8evt->Init(py8tree);
    std::cout << "CaloXPrimaryGeneratorAction: nentries=" << py8tree->GetEntriesFast() << std::endl;
  }
  else
  {
    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);

    // default particle kinematic
    //
    auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fParticleGun->SetParticleEnergy(50. * MeV);
  }
}

CaloXPrimaryGeneratorAction::~CaloXPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void CaloXPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  // cout<<"CaloXPrimaryGeneratorAction::GeneratePrimaries is called..."<<endl;
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //
  worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box *worldBox = nullptr;
  if (worldLV)
  {
    worldBox = dynamic_cast<G4Box *>(worldLV->GetSolid());
  }

  if (worldBox)
  {
    worldZHalfLength = worldBox->GetZHalfLength();
  }
  else
  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("CaloXPrimaryGeneratorAction::GeneratePrimaries()",
                "MyCode0002", JustWarning, msg);
  }

  // The size of Calorimeter volume

  double calorimeterZHalfLength = 0.;
  auto calorimeterLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Calorimeter");

  // Check that the world volume has box shape
  G4Box *calorimeterBox = nullptr;
  if (calorimeterLV)
  {
    calorimeterBox = dynamic_cast<G4Box *>(calorimeterLV->GetSolid());
  }

  if (calorimeterBox)
  {
    calorimeterZHalfLength = calorimeterBox->GetZHalfLength();
  }
  else
  {
    G4ExceptionDescription msg;
    msg << "Calorimeter volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("CaloXPrimaryGeneratorAction::GeneratePrimaries()",
                "MyCode0002", JustWarning, msg);
  }

  if (CaloXPythiaON == 1)
  {
    getPy8Event(anEvent);
  }
  else
  {
    float x = ((hh->getParamF("gun_x_max") - hh->getParamF("gun_x_min")) * G4UniformRand() + hh->getParamF("gun_x_min")) * cm;
    float y = ((hh->getParamF("gun_y_max") - hh->getParamF("gun_y_min")) * G4UniformRand() + hh->getParamF("gun_y_min")) * cm;
    float z = ((hh->getParamF("gun_z_max") - hh->getParamF("gun_z_min")) * G4UniformRand() + hh->getParamF("gun_z_min")) * cm;

    float en = ((hh->getParamF("gun_energy_max") - hh->getParamF("gun_energy_min")) * G4UniformRand() + hh->getParamF("gun_energy_min")) * GeV;
    string ptype = hh->getParamS("gun_particle");
    float px = ((hh->getParamF("pMomentum_x_max") - hh->getParamF("pMomentum_x_min")) * G4UniformRand() + hh->getParamF("pMomentum_x_min"));
    float py = ((hh->getParamF("pMomentum_y_max") - hh->getParamF("pMomentum_y_min")) * G4UniformRand() + hh->getParamF("pMomentum_y_min"));
    float pz = ((hh->getParamF("pMomentum_z_max") - hh->getParamF("pMomentum_z_min")) * G4UniformRand() + hh->getParamF("pMomentum_z_min"));

    // Re-use existing particle gun; update kinematics for this event
    auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(ptype);
    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
    fParticleGun->SetParticleEnergy(en);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
    fParticleGun->GeneratePrimaryVertex(anEvent);

    int pdgid = particleDefinition->GetPDGEncoding();
    hh->saveBeamXYZEPxPyPz(ptype, pdgid, x, y, z, en, px, py, pz);
  }
}

void CaloXPrimaryGeneratorAction::getParamFromEnvVars()
{
  CaloXPythiaON = hh->getParamB("CaloXPythiaON", false, 0); // default is 0, no Pythia8.

  CaloXPythiaXmin = hh->getParamD("CaloXPythiaXmin", false, 0.0);
  CaloXPythiaXmax = hh->getParamD("CaloXPythiaXmax", false, 0.0);
  CaloXPythiaYmin = hh->getParamD("CaloXPythiaYmin", false, 0.0);
  CaloXPythiaYmax = hh->getParamD("CaloXPythiaYmax", false, 0.0);
  CaloXPythiaZmin = hh->getParamD("CaloXPythiaZmin", false, -99999.0);
  CaloXPythiaZmax = hh->getParamD("CaloXPythiaZmax", false, CaloXPythiaZmin); // use the World boundary.

  CaloXPythiaSkip = hh->getParamI("CaloXPythiaSkip", false, 0);   // skip this many events.
  CaloXPythiaPrint = hh->getParamI("CaloXPythiaPrint", false, 5); // print this many events.

  CaloXPythiaFile = hh->getParamS("CaloXPythiaFile", false, "test.root");

  return;
}

// -----------------------------------------------------------------------------
void CaloXPrimaryGeneratorAction::getPy8Event(G4Event *anEvent)
{

  py8evt->GetEntry(py8eventNumber);
  if (py8eventCounter < CaloXPythiaPrint)
    printPy8Event();

  py8eventCounter++;
  py8eventNumber++;
  // std::cout<<"CaloXPrimaryGeneratorAction::getPy8Event  after py8evt->GetEntry="<<std::endl;
  // std::cout<<"py8evt->pid->size()  "<<py8evt->pid->size()<<std::endl;
  // std::cout<<"    pid=py8evt->pid->at(i) ="<<py8evt->pid->at(0)<<std::endl;

  double x = CaloXPythiaXmin + (CaloXPythiaXmax - CaloXPythiaXmin) * CLHEP::RandFlat::shoot();
  double y = CaloXPythiaYmin + (CaloXPythiaYmax - CaloXPythiaYmin) * CLHEP::RandFlat::shoot();
  double z = CaloXPythiaZmin + (CaloXPythiaZmax - CaloXPythiaZmin) * CLHEP::RandFlat::shoot();
  if (z < -worldZHalfLength)
    z = -worldZHalfLength + 0.0001; // limit to the World volume.
  double t = 0.0;

  // std::cout<<"CaloXPrimaryGeneratorAction::getPy8Event  x="<<x
  //<<"   CaloXPythiaXmin "<<CaloXPythiaXmin<<"  max "<<CaloXPythiaXmax<<std::endl;

  G4PrimaryVertex *vertex = new G4PrimaryVertex(G4ThreeVector(x, y, z), t);

  for (int i = 0; i < py8evt->particle_ID->size(); i++)
  {
    int pid = py8evt->particle_ID->at(i);
    double px = py8evt->particle_px->at(i) * GeV;
    double py = py8evt->particle_py->at(i) * GeV;
    double pz = py8evt->particle_pz->at(i) * GeV;
    G4ParticleDefinition *particle_definition = particleTable->FindParticle(pid);
    G4PrimaryParticle *particle = new G4PrimaryParticle(particle_definition, px, py, pz);
    vertex->SetPrimary(particle);
  } // end of loop over particles...

  anEvent->AddPrimaryVertex(vertex);

  // std::cout<<"anEvent->GetPrimaryVertex(0)->GetNumberOfParticle(): "<<anEvent->GetPrimaryVertex(0)->GetNumberOfParticle()<<std::endl;
  // std::cout<<"anEvent->GetPrimaryVertex(0)->GetZ0(): "<<anEvent->GetPrimaryVertex(0)->GetZ0()<<std::endl;
  // G4PrimaryParticle* primary = anEvent->GetPrimaryVertex(0)->GetPrimary(0);
  // std::cout<<"CaloXPrimaryGeneratorAction::getPy8Event:"<<primary->GetMomentum()<<std::endl;
}

void CaloXPrimaryGeneratorAction::printPy8Event()
{
  // Header.
  cout << "\n --------  PY8CaloX Event Listing ----------"
       << "-------------------------------------------------\n \n    no    "
       << "ID      p_x        p_y        p_z         e \n";
  cout << endl;

  // Precision. At high energy switch to scientific format for momenta.
  int precision = 3; // use 3 for now (sk)
  int prec = max(3, precision);

  for (int i = 0; i < py8evt->particle_E->size(); i++)
  {
    cout << setw(6) << std::right << i
         << setw(10) << std::right << py8evt->particle_ID->at(i)
         << setw(8 + prec) << py8evt->particle_px->at(i)
         << setw(8 + prec) << py8evt->particle_py->at(i)
         << setw(8 + prec) << py8evt->particle_pz->at(i)
         << setw(8 + prec) << py8evt->particle_E->at(i)
         << "\n";
  } // end of loop over particles...
}

void CaloXPrimaryGeneratorAction::FillHEPparticles(
    std::vector<int> *mHepPID, std::vector<int> *mHepStatus,
    std::vector<int> *mHepMother1, std::vector<int> *mHepMother2,
    std::vector<int> *mHepDaughter1, std::vector<int> *mHepDaughter2,
    std::vector<float> *mHepPx, std::vector<float> *mHepPy,
    std::vector<float> *mHepPz, std::vector<float> *mHepE,
    std::vector<float> *mHepMass)
{
  // NOTE: status, mother, daughter, and mass vectors are accepted in the signature
  // but not yet filled — Py8Jet does not currently store those branches.
  // Extend Py8Jet and add push_back calls here if that information is needed.
  (void)mHepStatus; (void)mHepMother1; (void)mHepMother2;
  (void)mHepDaughter1; (void)mHepDaughter2; (void)mHepMass;
  if (py8evt == NULL)
    return;

  int n = py8evt->particle_ID->size();
  for (int i = 0; i < n; i++)
  {
    mHepPID->push_back(py8evt->particle_ID->at(i));
    mHepPx->push_back(py8evt->particle_px->at(i));
    mHepPy->push_back(py8evt->particle_py->at(i));
    mHepPz->push_back(py8evt->particle_pz->at(i));
    mHepE->push_back(py8evt->particle_E->at(i));
  } // end of loop over particles...
}
