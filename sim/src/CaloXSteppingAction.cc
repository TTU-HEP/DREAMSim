/// \file CaloXSteppingAction.cc
/// \brief Implementation of the CaloXSteppingAction class

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"

#include "CaloXSteppingAction.hh"
#include "CaloXDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4Triton.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"

#include "CaloXID.h"  // including CaloXID, CaloXHit, CaloXTree
#include "CaloXHit.h" // including CaloXID, CaloXHit, CaloXTree
#include "CaloXTree.h"
#include "CaloXPhotonInfo.h"

#include "G4LogicalVolumeStore.hh"
#include "G4Tubs.hh"
#include "G4MaterialPropertiesTable.hh"


CaloXSteppingAction::CaloXSteppingAction(CaloXEventAction *eventAction, CaloXTree *histo)
    : G4UserSteppingAction(),
      fEventAction(eventAction),
      hh(histo)
{
  // initialize SiPM PDE
  int sipmType = histo->getParamI("sipmType");
  initPDE(sipmType); // 1= J 6 mm 6.0V, 2= J 6 mm 2.5V

  // Read optical photon sampling rod/layer from mac parameters (defaults: 45, 40)
  opSampleRod   = histo->getParamI("opSampleRod",   false, 45);
  opSampleLayer = histo->getParamI("opSampleLayer", false, 40);
  std::cout << "Optical photon sampling: rod=" << opSampleRod
            << "  layer=" << opSampleLayer << std::endl;



  std::cout << "  " << std::endl;
  std::cout << "sipmType " << sipmType << " in CaloXSteppingAction::CaloXSteppingAction" << std::endl;
  for (int i = 200; i < 900; i = i + 10)
  {
    float lambda = float(i) + 0.5;
    float pde = getPDE(lambda);
    std::cout << "i=" << i << "  lambda= " << lambda << "   pde= " << pde << std::endl;
  }
}


CaloXSteppingAction::~CaloXSteppingAction()
{
}


void CaloXSteppingAction::UserSteppingAction(const G4Step *step)
{
  G4Track *track = step->GetTrack();
  G4int trackID = track->GetTrackID();

  const G4VProcess *process = step->GetPostStepPoint()->GetProcessDefinedStep();
  G4String processName = "Unknown";
  if (process)
  {
    processName = process->GetProcessName();
  }

  // Collect energy and track length step by step

  //   === begin of checking optical photon ===
  static G4ParticleDefinition *opticalphoton =
      G4OpticalPhoton::OpticalPhotonDefinition();

  const G4DynamicParticle *dynamicParticle = track->GetDynamicParticle();
  const G4ParticleDefinition *particleDef =
      dynamicParticle->GetParticleDefinition();

  if (particleDef == opticalphoton)
  {
    fillOPInfo(step, false);
  }
  //   === end of checking optical photon ===

  // get volume of the current step
  // auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // energy deposit
  // edep does not include energy transfrerred to 2ndaries.
  // https://geant4-forum.web.cern.ch/t/total-energy-and-total-energy-deposit/6936/8
  auto edep = step->GetTotalEnergyDeposit();
  auto edepNonIon = step->GetNonIonizingEnergyDeposit();
  double charge = track->GetDefinition()->GetPDGCharge();

  G4int pdgcode = dynamicParticle->GetPDGcode();
  // G4int absPdgCode=abs(pdgcode);
  G4ParticleDefinition *particle = dynamicParticle->GetDefinition();
  G4String particleName = particle->GetParticleName();
  G4double kinEnergy = dynamicParticle->GetKineticEnergy();
  G4double mass = particle->GetPDGMass();

  //   energy deposit in cell...
  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto depth = touchable->GetHistory()->GetDepth();
  auto thisPhysical = touchable->GetVolume(); // mother
  auto thisCopyNo = thisPhysical->GetCopyNo();
  auto thisName = thisPhysical->GetName();
  G4ThreeVector posA = step->GetPreStepPoint()->GetPosition();

  double birks = 1.0;
  vector<double> ncer;

  int caloType = 0;
  int fiberNumber = 0;
  int holeNumber = 0;
  int rodNumber = 0;
  int layerNumber = 0;
  // int holeReplicaNumber=0;
  // int rodReplicaNumber=0;
  // int  layerReplicaNumber=0;

  // goes outside the world volume
  // condition taken from https://github.com/Geant4/geant4/blob/v11.2.2/examples/advanced/lAr_calorimeter/src/FCALSteppingAction.cc#L185
  if (track->GetNextVolume() == 0)
  {
    // std::cout<<"Stepping Action:  track goes outside the world volume"<<std::endl;
    double eLeak = step->GetPostStepPoint()->GetKineticEnergy();
    if (particle == G4Positron::Positron())
    {
      // electron pair production
      // mass is created from the vacuum using the kinetic energy
      // therefore should include it
      eLeak += 2 * electron_mass_c2;
    }
    if (posA.z() / cm <= 100 / cm)
    {
      hh->accumulateEnergy(eLeak / GeV, -92);
    }
    if (posA.z() / cm > 100 / cm)
    {
      hh->accumulateEnergy(eLeak / GeV, -91);
    }
    hh->accumulateEnergy(eLeak / GeV, -99);
  }

  // Track invisible energy (nuclear binding energy losses etc.)
  double e_net_change = findInvisible(step, false);
  hh->accumulateEnergy(e_net_change / GeV, -90);

  if (thisName.compare(0, 5, "World") == 0)
  {
    // outside the volume
    caloType = -1;
  }

  if (thisName.compare(0, 3, "Rod") == 0)
  {
    caloType = 1;
    fiberNumber = -1;
    holeNumber = 0;
    rodNumber = touchable->GetCopyNumber(0);
    layerNumber = touchable->GetCopyNumber(1);
    // holeReplicaNumber=touchable->GetReplicaNumber(2);
    // rodReplicaNumber=touchable->GetReplicaNumber(3);
    // layerReplicaNumber=touchable->GetReplicaNumber(4);
  }
  if (thisName.compare(0, 18, "fiberCoreScintPhys") == 0)
  {
    caloType = 2;
    birks = getBirk(step);
  }
  if (thisName == "fiberCoreCherePlasticPhys")
  {
    caloType = 3;
    ncer = UserCerenkov(step); // cerenkov photons (plastic fiber);
  }
  if (thisName == "fiberCoreChereQuartzPhys")
  {
    caloType = 4;
    ncer = UserCerenkov(step); // cerenkov photons (quartz fiber);
  }

  if (caloType == 2 || caloType == 3 || caloType == 4)
  {
    // todo: there might be better ways to get these information.
    // 0 is core number,
    // 1 is fiber/clad nummber (since core is in clad)
    // 2 is hole, 3 is rod, 4 is layer
    fiberNumber = touchable->GetCopyNumber(1);
    holeNumber = touchable->GetCopyNumber(2);
    rodNumber = touchable->GetCopyNumber(3);
    layerNumber = touchable->GetCopyNumber(4);
    // holeReplicaNumber=touchable->GetReplicaNumber(2);
    // rodReplicaNumber=touchable->GetReplicaNumber(3);
    // layerReplicaNumber=touchable->GetReplicaNumber(4);
  }

  hh->accumulateEnergy(edep / GeV, caloType);

  CaloXID caloid(caloType, fiberNumber, layerNumber, rodNumber, posA.z(), track->GetGlobalTime());

  CaloXHit aHit;
  aHit.caloid = caloid;
  aHit.x = posA.x() / cm; // in cm
  aHit.y = posA.y() / cm;
  aHit.z = posA.z() / cm;
  aHit.pid = pdgcode;
  aHit.fiberNumber = fiberNumber;
  aHit.rodNumber = rodNumber;
  aHit.layerNumber = layerNumber;
  // if (fabs(pdgcode) > 1e9)
  //{
  //   // https://indico.ph.tum.de/event/3955/sessions/752/attachments/2741/3099/Day3_Physics.pdf
  //   // larger than 1e9 is nuclei, ion, or excited state of atoms and nuclei
  //   // e.g., 1000020040 is alpha, 1000010020 is deuteron, 100ZZZAAAI* is
  //   // *ZZZ=proton number, AAA=nucleon number, I=excitation level
  //   std::cout << "pdgcode " << pdgcode << "  mass " << mass << "  kinEnergy " << kinEnergy << std::endl;
  // }
  aHit.calotype = caloType;
  aHit.trackid = track->GetTrackID();
  aHit.globaltime = track->GetGlobalTime() / ns;
  aHit.localtime = track->GetLocalTime() / ns;
  aHit.steplength = track->GetTrackLength() / cm;
  aHit.edep = edep / GeV;             //  in GeV
  aHit.edepNonIon = edepNonIon / GeV; // in GeV
  aHit.edepbirk = edep * birks / GeV;
  aHit.process = processName;
  if (ncer.size() > 0)
  {
    aHit.ncer = ncer[0];
    aHit.ncercap = ncer[3]; // including SiPM pde and capturing efficiency
  }

  // aHit.print();

  hh->accumulateHits(aHit);

  // hh->histo1D["edepX"]->Fill(aHit.x/10.0,edep);
  // hh->histo1D["edepY"]->Fill(aHit.y/10.0,edep);
  // hh->histo1D["cerX"]->Fill(aHit.x/10.0,ncer[0]);
  // hh->histo1D["cerY"]->Fill(aHit.y/10.0,ncer[0]);

  // fEventAction->AccumulateCaloXHits(aHit);
  // fEventAction->StepAnalysisSensor(step,ncer);   // analysis for the sensor volume.

  // fEventAction->StepAnalysis(step,ncer[0],ncer[1]);   // in original sim. Moved to above.
} // end of CaloXSteppingAction::UserSteppingAction.


vector<double> CaloXSteppingAction::UserCerenkov(const G4Step *step)
{
  double n_scint = 0;
  double n_cer = 0;
  double n_cerhad = 0;
  double nCERtotal = 0;
  double nCERlocal = 0;
  double nCERlocalElec = 0;
  double nCERlocalCap = 0;
  double nCERlocalElecCap = 0;

  //  Code from examples/extended/optical/OpNovice2/src/SteppingAction.cc
  static G4ParticleDefinition *opticalphoton =
      G4OpticalPhoton::OpticalPhotonDefinition();

  G4Track *track = step->GetTrack();
  G4StepPoint *endPoint = step->GetPostStepPoint();
  G4StepPoint *startPoint = step->GetPreStepPoint();

  const G4DynamicParticle *theParticle = track->GetDynamicParticle();
  const G4ParticleDefinition *particleDef =
      theParticle->GetParticleDefinition();
  G4int pdgcode = abs(theParticle->GetPDGcode());

  if (particleDef != opticalphoton)
  { // particle != opticalphoton
    // print how many Cerenkov and scint photons produced this step
    // this demonstrates use of GetNumPhotons()
    auto proc_man =
        track->GetDynamicParticle()->GetParticleDefinition()->GetProcessManager();
    G4ProcessVector *proc_vec = proc_man->GetPostStepProcessVector(typeDoIt);
    G4int n_proc = proc_vec->entries();

    // G4int n_scint = 0;
    // G4int n_cer   = 0;
    for (G4int i = 0; i < n_proc; ++i)
    {
      G4String proc_name = (*proc_vec)[i]->GetProcessName();
      if (proc_name.compare("Cerenkov") == 0)
      {
        auto cer = (G4Cerenkov *)(*proc_vec)[i];
        n_cer = cer->GetNumPhotons();
        if (pdgcode == 11)
        {
          n_cerhad = n_cer;
        }
      }
      else if (proc_name.compare("Scintillation") == 0)
      {
        auto scint = (G4Scintillation *)(*proc_vec)[i];
        n_scint = scint->GetNumPhotons();
      }
    }
    int fVerbose = -1; //  set a value here for now...
    if (fVerbose > 0)
    {
      if (n_cer > 0 || n_scint > 0)
      {
        std::cout << "In this step, " << n_cer << " Cerenkov and " << n_scint
                  << " scintillation photons were produced." << std::endl;
      }
    }
  }

  // loop over secondaries, create statistics
  const std::vector<const G4Track *> *secondaries =
      step->GetSecondaryInCurrentStep();

  for (auto sec : *secondaries)
  {
    if (sec->GetDynamicParticle()->GetParticleDefinition() == opticalphoton)
    {
      G4String creator_process = sec->GetCreatorProcess()->GetProcessName();
      if (creator_process.compare("Cerenkov") == 0)
      {
        G4double en = sec->GetKineticEnergy();
        double wavelength = 1239.8 * eV / en;
        G4ThreeVector pvec = sec->GetMomentumDirection();

        int capture = 0;
        if (abs(pvec.theta()) < 0.336)
          capture = 1; // NA=sin(theta)=0.33

        nCERtotal = nCERtotal + 1;

        double pde = getPDE(wavelength);
        nCERlocal = nCERlocal + pde;

        if (capture == 1)
        {
          nCERlocalCap = nCERlocalCap + pde;
        }

        if (pdgcode == 11)
        {
          nCERlocalElec = nCERlocalElec + pde;
          if (capture == 1)
          {
            nCERlocalElecCap = nCERlocalElecCap + pde;
          }
        }
        // cout<<"cerenkov phton  en="<<en<<endl;
        // run->AddCerenkovEnergy(en);
        // run->AddCerenkov();
        // analysisMan->FillH1(1, en / eV);
      }
      else if (creator_process.compare("Scintillation") == 0)
      {
        G4double en = sec->GetKineticEnergy();
        // run->AddScintillationEnergy(en);
        // run->AddScintillation();
        // analysisMan->FillH1(2, en / eV);

        // G4double time = sec->GetGlobalTime();
        // analysisMan->FillH1(3, time / ns);
      }
    }
  } //  end of for(auto sec : *secondaries)

  // double NCER=double(n_cer)/10000.0;
  vector<double> NCER;
  NCER.push_back(double(nCERtotal));
  NCER.push_back(double(nCERlocal));
  NCER.push_back(double(nCERlocalElec));
  NCER.push_back(double(nCERlocalCap));
  NCER.push_back(double(nCERlocalElecCap));
  return NCER;
}

// ========================================================================================
double CaloXSteppingAction::getBirk(const G4Step *step)
{
  double weight = 1.0;

  double edep = step->GetTotalEnergyDeposit();
  double steplength = step->GetStepLength() / CLHEP::cm; //  convert from mm to cm.
  double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
  double density = step->GetPreStepPoint()->GetMaterial()->GetDensity() / (CLHEP::g / CLHEP::cm3);
  string materialName = step->GetPreStepPoint()->GetMaterial()->GetName();
  // std::cout<<"CaloXEventAction::getBirk:   name="<<materialName<<std::endl;

  /*
       if(materialName.compare(0,8,"G4_PbWO4")==0) {
          weight=getBirkL3(edep,steplength,charge,density);
          // std::cout<<"CaloXEventAction::getBirk:   name="<<materialName<<"   weight="<<weight<<std::endl;
       }

       if(materialName.compare(0,14,"H_Scintillator")==0 || materialName.compare(0,26,"G4_PLASTIC_SC_VINYLTOLUENE")==0) {
          weight=getBirkHC(edep,steplength,charge,density);
          // std::cout<<"CaloXEventAction::getBirk:   name="<<materialName<<"   weight="<<weight<<std::endl;
       }
  */

  if (materialName.compare(0, 11, "Polystyrene") == 0)
  {
    weight = getBirkHC(edep, steplength, charge, density);
    // std::cout<<"CaloXEventAction::getBirk:   name="<<materialName<<"   weight="<<weight<<std::endl;
  }

  return weight;
}

double CaloXSteppingAction::getBirkHC(double dEStep, double step, double charge, double density)
{
  double weight = 1.;
  if (charge != 0. && step > 0.)
  {
    double birkC1HC_ = 0.0052;
    double birkC2HC_ = 0.142;
    double birkC3HC_ = 1.75;
    double dedx = dEStep / step;
    double rkb = birkC1HC_ / density;
    double c = birkC2HC_ * rkb * rkb;
    if (std::abs(charge) >= 2.)
      rkb /= birkC3HC_;
    weight = 1. / (1. + rkb * dedx + c * dedx * dedx);
  }
  return weight;
}

double CaloXSteppingAction::getBirkL3(double dEStep, double step, double charge, double density)
{
  double weight = 1.;
  if (charge != 0. && step > 0.)
  {
    double birkC1EC_ = 0.03333;
    double birkSlopeEC_ = 0.253694;
    double birkCutEC_ = 0.1;
    double dedx = dEStep / step;
    double rkb = birkC1EC_ / density;
    if (dedx > 0)
    {
      weight = 1. - birkSlopeEC_ * log(rkb * dedx);
      if (weight < birkCutEC_)
        weight = birkCutEC_;
      else if (weight > 1.)
        weight = 1.;
    }
  }
  return weight;
}

void CaloXSteppingAction::initPDE(int sipmType)
{
  //  sipmType  1= J-6mm-6.0V,  2=J-6mm-2.5V
  //
  //  700 bins of wavelength from 200 to 899 nm
  pdeLambdaMin = 200;
  pdeLambdaMax = 899;
  //
  //
  //
  std::vector<float> pde_J_6mm_6p0 = {
      4.75, 4.75, 4.75, 4.76, 4.77, 4.78, 4.79, 4.81, 4.82, 4.84,
      4.86, 4.87, 4.89, 4.90, 4.91, 4.92, 4.92, 4.92, 4.92, 4.92,
      4.91, 4.91, 4.91, 4.91, 4.91, 4.91, 4.92, 4.92, 4.92, 4.92,
      4.92, 4.92, 4.91, 4.91, 4.90, 4.91, 4.91, 4.93, 4.95, 4.99,
      5.03, 5.08, 5.14, 5.20, 5.27, 5.34, 5.41, 5.47, 5.51, 5.54,
      5.56, 5.56, 5.57, 5.57, 5.58, 5.59, 5.62, 5.67, 5.73, 5.81,
      5.90, 6.00, 6.12, 6.25, 6.40, 6.56, 6.73, 6.94, 7.17, 7.44,
      7.75, 8.10, 8.48, 8.87, 9.26, 9.65, 10.09, 10.60, 11.23, 11.93,
      12.65, 13.28, 13.83, 14.29, 14.71, 15.11, 15.53, 15.99, 16.53, 17.16,
      17.86, 18.55, 19.20, 19.87, 20.60, 21.43, 22.36, 23.26, 23.94, 24.38,
      24.67, 24.89, 25.13, 25.44, 25.80, 26.20, 26.61, 27.02, 27.40, 27.74,
      28.06, 28.36, 28.64, 28.92, 29.20, 29.50, 29.81, 30.13, 30.45, 30.78,
      31.10, 31.41, 31.71, 32.01, 32.29, 32.56, 32.82, 33.07, 33.31, 33.54,
      33.77, 34.01, 34.25, 34.50, 34.75, 35.00, 35.23, 35.44, 35.61, 35.75,
      35.85, 35.93, 36.00, 36.06, 36.13, 36.23, 36.34, 36.46, 36.59, 36.72,
      36.86, 36.98, 37.10, 37.21, 37.30, 37.39, 37.48, 37.57, 37.66, 37.75,
      37.85, 37.96, 38.08, 38.21, 38.34, 38.49, 38.65, 38.82, 39.01, 39.22,
      39.45, 39.71, 40.00, 40.33, 40.69, 41.07, 41.46, 41.83, 42.15, 42.42,
      42.65, 42.86, 43.06, 43.25, 43.46, 43.69, 43.93, 44.18, 44.43, 44.69,
      44.95, 45.21, 45.45, 45.68, 45.88, 46.05, 46.20, 46.33, 46.45, 46.58,
      46.71, 46.85, 47.02, 47.21, 47.41, 47.62, 47.83, 48.03, 48.22, 48.40,
      48.56, 48.71, 48.86, 48.99, 49.11, 49.23, 49.33, 49.43, 49.52, 49.60,
      49.67, 49.73, 49.78, 49.81, 49.84, 49.85, 49.85, 49.85, 49.83, 49.81,
      49.77, 49.73, 49.69, 49.64, 49.58, 49.52, 49.46, 49.39, 49.31, 49.23,
      49.14, 49.04, 48.93, 48.82, 48.70, 48.57, 48.45, 48.31, 48.18, 48.04,
      47.90, 47.76, 47.61, 47.45, 47.29, 47.12, 46.95, 46.77, 46.58, 46.39,
      46.20, 46.00, 45.80, 45.60, 45.39, 45.19, 45.00, 44.80, 44.62, 44.44,
      44.27, 44.10, 43.94, 43.77, 43.61, 43.45, 43.29, 43.12, 42.94, 42.76,
      42.57, 42.37, 42.16, 41.93, 41.70, 41.44, 41.18, 40.92, 40.67, 40.43,
      40.22, 40.03, 39.86, 39.68, 39.50, 39.30, 39.07, 38.80, 38.52, 38.23,
      37.95, 37.70, 37.48, 37.27, 37.07, 36.85, 36.61, 36.33, 36.03, 35.71,
      35.40, 35.10, 34.84, 34.61, 34.42, 34.24, 34.09, 33.93, 33.78, 33.62,
      33.44, 33.24, 33.02, 32.79, 32.57, 32.35, 32.15, 31.96, 31.79, 31.62,
      31.44, 31.25, 31.05, 30.81, 30.56, 30.29, 30.03, 29.78, 29.55, 29.35,
      29.19, 29.05, 28.93, 28.81, 28.70, 28.57, 28.43, 28.28, 28.11, 27.94,
      27.75, 27.57, 27.38, 27.19, 27.01, 26.83, 26.65, 26.47, 26.30, 26.13,
      25.96, 25.80, 25.64, 25.48, 25.32, 25.15, 24.98, 24.80, 24.63, 24.47,
      24.31, 24.17, 24.03, 23.90, 23.77, 23.64, 23.51, 23.37, 23.22, 23.08,
      22.93, 22.79, 22.65, 22.51, 22.39, 22.28, 22.18, 22.10, 22.01, 21.93,
      21.84, 21.75, 21.64, 21.51, 21.37, 21.21, 21.04, 20.86, 20.68, 20.51,
      20.34, 20.18, 20.03, 19.90, 19.78, 19.66, 19.56, 19.47, 19.38, 19.29,
      19.21, 19.12, 19.03, 18.94, 18.84, 18.74, 18.64, 18.53, 18.41, 18.30,
      18.18, 18.06, 17.94, 17.82, 17.71, 17.59, 17.48, 17.36, 17.25, 17.14,
      17.03, 16.93, 16.83, 16.73, 16.63, 16.53, 16.43, 16.34, 16.24, 16.14,
      16.04, 15.94, 15.84, 15.74, 15.63, 15.53, 15.43, 15.32, 15.22, 15.12,
      15.02, 14.92, 14.83, 14.74, 14.65, 14.56, 14.48, 14.40, 14.31, 14.24,
      14.16, 14.08, 14.00, 13.93, 13.85, 13.77, 13.69, 13.61, 13.53, 13.44,
      13.36, 13.27, 13.19, 13.10, 13.01, 12.92, 12.83, 12.75, 12.66, 12.57,
      12.49, 12.40, 12.32, 12.23, 12.15, 12.07, 11.99, 11.91, 11.82, 11.74,
      11.66, 11.57, 11.49, 11.40, 11.32, 11.23, 11.14, 11.05, 10.96, 10.87,
      10.77, 10.68, 10.58, 10.48, 10.38, 10.28, 10.18, 10.07, 9.97, 9.86,
      9.76, 9.65, 9.55, 9.46, 9.36, 9.28, 9.19, 9.11, 9.04, 8.98,
      8.91, 8.86, 8.80, 8.75, 8.71, 8.66, 8.62, 8.58, 8.54, 8.51,
      8.47, 8.43, 8.39, 8.36, 8.32, 8.27, 8.23, 8.18, 8.13, 8.08,
      8.02, 7.95, 7.89, 7.82, 7.75, 7.68, 7.60, 7.53, 7.45, 7.38,
      7.30, 7.23, 7.16, 7.09, 7.02, 6.96, 6.90, 6.84, 6.78, 6.72,
      6.66, 6.60, 6.55, 6.49, 6.44, 6.39, 6.33, 6.28, 6.23, 6.17,
      6.12, 6.07, 6.02, 5.97, 5.92, 5.88, 5.83, 5.79, 5.75, 5.71,
      5.68, 5.65, 5.62, 5.59, 5.55, 5.52, 5.49, 5.46, 5.42, 5.39,
      5.35, 5.30, 5.26, 5.21, 5.16, 5.11, 5.06, 5.00, 4.95, 4.91,
      4.86, 4.82, 4.78, 4.74, 4.71, 4.68, 4.65, 4.63, 4.60, 4.57,
      4.55, 4.52, 4.49, 4.46, 4.43, 4.39, 4.35, 4.30, 4.26, 4.22,
      4.17, 4.13, 4.09, 4.05, 4.02, 4.00, 3.97, 3.95, 3.94, 3.93,
      3.91, 3.90, 3.89, 3.87, 3.85, 3.82, 3.79, 3.76, 3.73, 3.69,
      3.65, 3.61, 3.57, 3.53, 3.49, 3.45, 3.41, 3.38, 3.34, 3.31,
      3.28, 3.26, 3.23, 3.20, 3.17, 3.15, 3.12, 3.09, 3.06, 3.03,
      3.00, 2.97, 2.94, 2.90, 2.87, 2.83, 2.80, 2.76, 2.72, 2.68,
      2.64, 2.60, 2.55, 2.51, 2.47, 2.42, 2.38, 2.34, 2.30, 2.26,
      2.22, 2.19, 2.16, 2.13, 2.10, 2.07, 2.05, 2.03, 2.01, 1.99,
      1.98, 1.97, 1.96, 1.96, 1.96, 1.97, 1.97, 1.98, 2.00, 2.02}; // end of pde_J_6mm_6p0=

  for (float &element : pde_J_6mm_6p0)
  {
    element *= 0.01;
  }

  //
  //
  //
  std::vector<float> pde_J_6mm_2p5 = {
      2.07, 2.14, 2.21, 2.29, 2.38, 2.46, 2.56, 2.65, 2.74, 2.84,
      2.93, 3.01, 3.09, 3.16, 3.22, 3.27, 3.31, 3.33, 3.35, 3.35,
      3.35, 3.34, 3.33, 3.31, 3.30, 3.28, 3.27, 3.26, 3.25, 3.24,
      3.24, 3.24, 3.24, 3.24, 3.24, 3.24, 3.23, 3.22, 3.21, 3.22,
      3.24, 3.28, 3.35, 3.44, 3.53, 3.63, 3.73, 3.82, 3.89, 3.93,
      3.94, 3.94, 3.93, 3.91, 3.91, 3.92, 3.95, 4.01, 4.09, 4.17,
      4.24, 4.31, 4.38, 4.46, 4.55, 4.67, 4.83, 5.03, 5.28, 5.59,
      5.93, 6.28, 6.64, 6.96, 7.24, 7.47, 7.68, 7.90, 8.14, 8.43,
      8.80, 9.26, 9.85, 10.51, 11.12, 11.60, 11.98, 12.37, 12.85, 13.40,
      13.94, 14.42, 14.82, 15.21, 15.67, 16.27, 17.03, 17.79, 18.40, 18.83,
      19.15, 19.42, 19.70, 20.06, 20.47, 20.91, 21.32, 21.67, 21.99, 22.27,
      22.55, 22.83, 23.14, 23.44, 23.75, 24.04, 24.30, 24.52, 24.72, 24.90,
      25.07, 25.24, 25.40, 25.58, 25.77, 25.99, 26.23, 26.50, 26.77, 27.03,
      27.28, 27.49, 27.67, 27.82, 27.96, 28.07, 28.18, 28.28, 28.39, 28.50,
      28.61, 28.72, 28.83, 28.95, 29.06, 29.17, 29.28, 29.38, 29.47, 29.56,
      29.65, 29.74, 29.82, 29.91, 30.00, 30.09, 30.20, 30.31, 30.43, 30.56,
      30.68, 30.81, 30.93, 31.04, 31.15, 31.24, 31.32, 31.40, 31.50, 31.61,
      31.75, 31.94, 32.18, 32.46, 32.78, 33.11, 33.43, 33.72, 33.97, 34.19,
      34.39, 34.58, 34.76, 34.96, 35.18, 35.40, 35.63, 35.85, 36.04, 36.19,
      36.30, 36.39, 36.44, 36.49, 36.52, 36.55, 36.59, 36.64, 36.70, 36.79,
      36.90, 37.01, 37.14, 37.26, 37.38, 37.50, 37.60, 37.70, 37.78, 37.85,
      37.91, 37.96, 38.00, 38.03, 38.06, 38.08, 38.09, 38.10, 38.09, 38.09,
      38.07, 38.05, 38.01, 37.97, 37.92, 37.86, 37.79, 37.71, 37.62, 37.52,
      37.42, 37.31, 37.20, 37.09, 36.97, 36.85, 36.73, 36.61, 36.49, 36.37,
      36.24, 36.11, 35.99, 35.86, 35.73, 35.59, 35.46, 35.32, 35.18, 35.03,
      34.88, 34.72, 34.57, 34.41, 34.25, 34.08, 33.92, 33.76, 33.60, 33.44,
      33.28, 33.12, 32.97, 32.81, 32.65, 32.48, 32.30, 32.11, 31.92, 31.70,
      31.48, 31.25, 31.02, 30.79, 30.57, 30.37, 30.18, 30.01, 29.86, 29.71,
      29.57, 29.44, 29.30, 29.17, 29.03, 28.89, 28.74, 28.58, 28.40, 28.21,
      28.00, 27.77, 27.54, 27.30, 27.06, 26.83, 26.60, 26.39, 26.19, 25.99,
      25.80, 25.62, 25.44, 25.27, 25.11, 24.95, 24.79, 24.63, 24.47, 24.31,
      24.14, 23.97, 23.81, 23.64, 23.48, 23.32, 23.17, 23.03, 22.90, 22.78,
      22.66, 22.55, 22.44, 22.32, 22.20, 22.07, 21.94, 21.80, 21.66, 21.52,
      21.38, 21.24, 21.12, 20.99, 20.85, 20.72, 20.57, 20.41, 20.23, 20.04,
      19.84, 19.64, 19.45, 19.28, 19.12, 18.99, 18.89, 18.80, 18.71, 18.64,
      18.57, 18.49, 18.40, 18.30, 18.20, 18.09, 17.97, 17.86, 17.75, 17.64,
      17.53, 17.43, 17.34, 17.25, 17.16, 17.07, 16.97, 16.86, 16.75, 16.63,
      16.51, 16.39, 16.27, 16.15, 16.04, 15.94, 15.84, 15.75, 15.66, 15.58,
      15.50, 15.42, 15.34, 15.25, 15.17, 15.08, 14.99, 14.90, 14.81, 14.71,
      14.61, 14.51, 14.41, 14.30, 14.20, 14.09, 13.99, 13.89, 13.79, 13.69,
      13.59, 13.50, 13.41, 13.33, 13.24, 13.16, 13.08, 13.00, 12.92, 12.84,
      12.76, 12.69, 12.61, 12.53, 12.45, 12.38, 12.30, 12.22, 12.14, 12.06,
      11.99, 11.91, 11.84, 11.76, 11.69, 11.63, 11.56, 11.50, 11.44, 11.38,
      11.32, 11.26, 11.20, 11.13, 11.06, 10.99, 10.91, 10.82, 10.74, 10.65,
      10.56, 10.47, 10.38, 10.30, 10.22, 10.14, 10.07, 10.00, 9.93, 9.88,
      9.82, 9.78, 9.73, 9.68, 9.64, 9.60, 9.55, 9.51, 9.46, 9.41,
      9.35, 9.29, 9.23, 9.17, 9.10, 9.04, 8.97, 8.91, 8.84, 8.78,
      8.72, 8.66, 8.60, 8.55, 8.49, 8.44, 8.39, 8.34, 8.29, 8.24,
      8.19, 8.14, 8.10, 8.05, 8.00, 7.95, 7.91, 7.86, 7.81, 7.76,
      7.71, 7.66, 7.61, 7.56, 7.50, 7.45, 7.39, 7.33, 7.27, 7.21,
      7.15, 7.09, 7.03, 6.97, 6.91, 6.85, 6.78, 6.72, 6.66, 6.60,
      6.54, 6.48, 6.43, 6.37, 6.32, 6.26, 6.21, 6.16, 6.11, 6.06,
      6.01, 5.96, 5.91, 5.86, 5.81, 5.76, 5.71, 5.66, 5.62, 5.58,
      5.55, 5.52, 5.49, 5.47, 5.45, 5.44, 5.44, 5.43, 5.42, 5.42,
      5.41, 5.39, 5.37, 5.33, 5.29, 5.24, 5.18, 5.11, 5.03, 4.95,
      4.87, 4.78, 4.69, 4.61, 4.53, 4.45, 4.37, 4.30, 4.24, 4.19,
      4.15, 4.12, 4.09, 4.06, 4.05, 4.03, 4.01, 4.00, 3.98, 3.96,
      3.94, 3.91, 3.88, 3.85, 3.81, 3.78, 3.74, 3.69, 3.65, 3.61,
      3.57, 3.52, 3.48, 3.44, 3.41, 3.37, 3.34, 3.31, 3.29, 3.26,
      3.24, 3.22, 3.21, 3.19, 3.18, 3.16, 3.15, 3.13, 3.12, 3.11,
      3.09, 3.08, 3.06, 3.05, 3.03, 3.01, 2.99, 2.97, 2.94, 2.92,
      2.89, 2.87, 2.84, 2.82, 2.79, 2.76, 2.73, 2.71, 2.68, 2.65,
      2.63, 2.60, 2.57, 2.55, 2.52, 2.49, 2.46, 2.43, 2.40, 2.38,
      2.35, 2.31, 2.28, 2.25, 2.22, 2.19, 2.15, 2.12, 2.08, 2.05,
      2.01, 1.98, 1.94, 1.91, 1.88, 1.85, 1.82, 1.79, 1.76, 1.74,
      1.72, 1.70, 1.68, 1.66, 1.65, 1.64, 1.64, 1.63, 1.63, 1.63,
      1.63, 1.63, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64,
      1.64, 1.63, 1.63, 1.63, 1.63, 1.63, 1.63, 1.63, 1.63, 1.63,
      1.64, 1.64, 1.63, 1.63, 1.63, 1.63, 1.63, 1.62, 1.61, 1.61,
      1.60, 1.59, 1.57, 1.56, 1.54, 1.52, 1.50, 1.47, 1.45, 1.42}; // end of pde_J_6mm_2p5=

  for (float &element : pde_J_6mm_2p5)
  {
    element *= 0.01;
  }

  if (sipmType == 1)
    sipmPDE = pde_J_6mm_6p0;
  if (sipmType == 2)
    sipmPDE = pde_J_6mm_2p5;
}

float CaloXSteppingAction::getPDE(float lambda)
{
  float pde = 0.0;
  int k = int(lambda);
  if (k > pdeLambdaMin and k < pdeLambdaMax)
  {
    int j = k - pdeLambdaMin;
    // std::cout<<"lambda "<<lambda<<"  k "<<k<<"  j "<<j<<std::endl;
    pde = sipmPDE[j];
  }
  return pde;
}

double CaloXSteppingAction::findInvisible(const G4Step *step, bool verbose)
{
  // check energy conservation
  // return the energy change that is not GetTotalEnergyDeposit in edep in Geant
  double e_kinetic_lost = step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy();
  double edep = step->GetTotalEnergyDeposit();
  double edepNonIon = step->GetNonIonizingEnergyDeposit();

  double e_secondary = 0.0;
  double m_produced = 0.0;
  int nProton = 0;
  int nNeutron = 0;

  const std::vector<const G4Track *> *secondaries = step->GetSecondaryInCurrentStep();
  for (auto sec : *secondaries)
  {
    e_secondary += sec->GetKineticEnergy();
    m_produced += sec->GetDynamicParticle()->GetMass();
    if (sec->GetParticleDefinition() == G4Proton::Proton())
    {
      nProton++;
    }
    if (sec->GetParticleDefinition() == G4Neutron::Neutron())
    {
      nNeutron++;
    }
  }

  // deal with the decay process,
  // where the energy of the post-step should be included
  // triton seems to be a special case, with zero secondary particles after the decay step
  G4ParticleDefinition *particle = step->GetTrack()->GetDefinition();
  if (step->GetPostStepPoint()->GetProcessDefinedStep() &&
      step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().compare("Decay") == 0 &&
      particle != G4Triton::Triton())
  {
    // in Geant4 at the decay step, it seems the total energy (i.e., mass, with kinetic energy beign zero)
    // is not reset to 0
    e_secondary -= step->GetPostStepPoint()->GetTotalEnergy();
    // mass of the decay product is created from the vacuum using the kinetic energy
    e_secondary += m_produced;
  }
  double e_net_change = e_kinetic_lost - e_secondary - edep;

  // subtract the mass of the produced particles for protons and neutrons
  // not totally sure if this is the right way to do it
  double mProton = G4Proton::ProtonDefinition()->GetPDGMass();
  while (e_net_change > mProton && nProton > 0)
  {
    e_net_change -= mProton;
    nProton--;
  }
  double mNeutron = G4Neutron::NeutronDefinition()->GetPDGMass();
  while (e_net_change > mNeutron && nNeutron > 0)
  {
    e_net_change -= mNeutron;
    nNeutron--;
  }

  if (verbose && fabs(e_net_change) >= 100.0 * MeV)
  {
    std::cout << "Found energy conservation problem: " << e_net_change / MeV << " MeV. Dump information.." << std::endl;
    std::cout << "Step number: " << step->GetTrack()->GetCurrentStepNumber()
              << ", Particle name: " << step->GetTrack()->GetDefinition()->GetParticleName()
              << ", Number of secondaries: " << secondaries->size()
              << ", Energy total lost: MeV " << e_net_change / MeV
              << ", Energy to secondaries: " << e_secondary / MeV
              << ", Energy deposited: " << edep / MeV
              << ", Energy deposited non-ionizing: " << edepNonIon / MeV
              << ", Delta Kinetic energy: " << e_kinetic_lost / MeV
              << ", Delta Total energy: " << step->GetPreStepPoint()->GetTotalEnergy() - step->GetPostStepPoint()->GetTotalEnergy()
              << std::endl;
    std::cout << ", Pre-step position: (" << step->GetPreStepPoint()->GetPosition().x() << ", " << step->GetPreStepPoint()->GetPosition().y() << ", " << step->GetPreStepPoint()->GetPosition().z() << ")"
              << ", Post-step position: (" << step->GetPostStepPoint()->GetPosition().x() << ", " << step->GetPostStepPoint()->GetPosition().y() << ", " << step->GetPostStepPoint()->GetPosition().z() << ")"
              << ", Pre-step time: " << step->GetPreStepPoint()->GetGlobalTime()
              << ", Post-step time: " << step->GetPostStepPoint()->GetGlobalTime()
              << ", Pre-step momentum: (" << step->GetPreStepPoint()->GetMomentum().x() << ", " << step->GetPreStepPoint()->GetMomentum().y() << ", " << step->GetPreStepPoint()->GetMomentum().z() << ") with kinetic energy: " << step->GetPreStepPoint()->GetKineticEnergy() << " and total energy: " << step->GetPreStepPoint()->GetTotalEnergy()
              << ", Post-step momentum: (" << step->GetPostStepPoint()->GetMomentum().x() << ", " << step->GetPostStepPoint()->GetMomentum().y() << ", " << step->GetPostStepPoint()->GetMomentum().z() << ") with kinetic energy: " << step->GetPostStepPoint()->GetKineticEnergy() << " and total energy: " << step->GetPostStepPoint()->GetTotalEnergy()
              << std::endl;
    if (step->GetPreStepPoint()->GetProcessDefinedStep())
    {
      std::cout << "Process name: (prestep) " << step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
    }
    if (step->GetPostStepPoint()->GetProcessDefinedStep())
    {
      std::cout << "Process name: (poststep) " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
    }
    double e_kin_sum = 0;
    double e_tot_sum = 0;
    for (auto sec : *secondaries)
    {
      e_kin_sum += sec->GetKineticEnergy();
      e_tot_sum += sec->GetTotalEnergy();
      std::cout << "Secondary name: " << sec->GetParticleDefinition()->GetParticleName()
                << ", Kinetic energy: " << sec->GetKineticEnergy() / MeV
                << ", Position: (" << sec->GetPosition().x() << ", " << sec->GetPosition().y() << ", " << sec->GetPosition().z() << ")"
                << ", Momentum: (" << sec->GetMomentum().x() << ", " << sec->GetMomentum().y() << ", " << sec->GetMomentum().z() << ")"
                << ", Process name: " << sec->GetCreatorProcess()->GetProcessName()
                << std::endl;
    }
    std::cout << "Sum of kinetic energy of secondaries: " << e_kin_sum / MeV << " MeV"
              << ", Sum of total energy of secondaries: " << e_tot_sum / MeV << " MeV"
              << std::endl;
    std::cout << "Track status " << step->GetTrack()->GetTrackStatus() << std::endl;
    std::cout << std::endl;
  }
  return e_net_change;
}

// ---------------------------------------------------------------------------
// Helper: compute both meridional and skew TIR results and store into photon.
// Called once at photon birth for every core photon, regardless of rod.
// ---------------------------------------------------------------------------
static void computeAnalyticalTIR(CaloXPhotonInfo &photon,
                                  const G4ThreeVector &pdir,
                                  const G4ThreeVector &pos,
                                  double n_core, double n_clad,
                                  double R_core, double fiber_half_z,
                                  double global_time,
                                  double abs_len)
{
  const double abs_pz      = std::abs(pdir.z());
  const double z_exit      = (pdir.z() >= 0.0) ? fiber_half_z : -fiber_half_z;
  const double z_travel    = std::abs(z_exit - pos.z());
  const double tir_thresh  = n_clad / n_core;

  // Arrival time for both methods uses abs_pz, NOT cos(captureAngle_s).
  // The axial speed is v_z = (c/n_core)*|pz| for any ray — skew or meridional.
  // captureAngle_s is the incidence angle at the cylindrical wall, which is a
  // different angle; cos(captureAngle_s) != |pz| for skew rays.
  // The two arrival times differ only in which photons qualify (TIR condition),
  // not in the propagation formula itself.

  // Total path length through the core along the ray (longer than z_travel for non-axial rays)
  const double path_length = (abs_pz > 0.0) ? z_travel / abs_pz : DBL_MAX;
  // Beer-Lambert survival probability: P = exp(-path / abs_len)
  // Sample once and apply the same random decision to both methods so that
  // isCaptured_m and isCaptured_s are comparable on a photon-by-photon basis.
  const double survival_prob = std::exp(-path_length / abs_len);
  const bool   survives      = (CLHEP::RandFlat::shoot() < survival_prob);

  // Meridional: treats fiber as slab, ignores azimuthal angular momentum
  photon.captureAngle_m = std::acos(std::min(abs_pz, 1.0));
  photon.isCaptured_m   = (abs_pz >= tir_thresh);
  photon.isAttenuated_m = photon.isCaptured_m && !survives;
  if (photon.isCaptured_m && survives && abs_pz > 0.0)
    photon.analyticalArrivalTime_m =
        (global_time + z_travel * n_core / (CLHEP::c_light * abs_pz)) / ns;

  // Skew-ray: uses conserved z-angular momentum L_z = x*py - y*px.
  // Wall incidence angle: sin(theta_wall) = sqrt(pz^2 + Lz^2/R^2) — used for TIR only.
  // captureAngle_s stores alpha = acos(|pz|), same definition as captureAngle_m,
  // because alpha is the angle that governs axial speed and travel time for both methods.
  const double Lz         = pos.x() * pdir.y() - pos.y() * pdir.x();
  const double sin_theta_wall = std::min(std::sqrt(pdir.z()*pdir.z() + Lz*Lz / (R_core*R_core)), 1.0);
  photon.captureAngle_s = std::acos(std::min(abs_pz, 1.0)); // alpha, same as _m
  photon.isCaptured_s   = (sin_theta_wall >= tir_thresh);
  photon.isAttenuated_s = photon.isCaptured_s && !survives;
  if (photon.isCaptured_s && survives && abs_pz > 0.0)
    photon.analyticalArrivalTime_s =
        (global_time + z_travel * n_core / (CLHEP::c_light * abs_pz)) / ns;
}

void CaloXSteppingAction::initOptics()
{
  // Called lazily on the first step so that materials and volumes are guaranteed
  // to exist (CaloXDetectorConstruction::Construct() has already run).
  const G4double refEnergy = 3.0 * eV; // mid-visible; RINDEX tables are flat here

  auto readRIndex = [&](const char *matName) -> double {
    G4Material *mat = G4Material::GetMaterial(matName);
    if (!mat)
      G4Exception("CaloXSteppingAction::initOptics", "NoMaterial", FatalException,
                  (G4String("Material not found: ") + matName).c_str());
    G4MaterialPropertiesTable *mpt = mat->GetMaterialPropertiesTable();
    if (!mpt)
      G4Exception("CaloXSteppingAction::initOptics", "NoMPT", FatalException,
                  (G4String("No MaterialPropertiesTable for: ") + matName).c_str());
    return mpt->GetProperty("RINDEX")->Value(refEnergy);
  };

  fN_CoreS = readRIndex("Polystyrene");
  fN_CladS = readRIndex("PMMA_Clad");
  fN_CoreC = readRIndex("PMMA");
  fN_CladC = readRIndex("Fluorinated_Polymer");
  fN_CoreQ = readRIndex("Fused_Silica");
  fN_CladQ = readRIndex("Hard_Polymer");

  // Read bulk absorption lengths (ABSLENGTH) at the same reference energy
  auto readAbsLen = [&](const char *matName) -> double {
    G4Material *mat = G4Material::GetMaterial(matName);
    G4MaterialPropertiesTable *mpt = mat->GetMaterialPropertiesTable();
    G4MaterialPropertyVector *vec = mpt->GetProperty("ABSLENGTH");
    return vec ? vec->Value(refEnergy) : DBL_MAX; // no table → no attenuation
  };
  fAbsLen_CoreS = readAbsLen("Polystyrene");
  fAbsLen_CoreC = readAbsLen("PMMA");
  fAbsLen_CoreQ = readAbsLen("Fused_Silica");

  auto *coreVol = G4LogicalVolumeStore::GetInstance()->GetVolume("fiberCoreS");
  if (!coreVol)
    G4Exception("CaloXSteppingAction::initOptics", "NoVolume", FatalException,
                "Logical volume fiberCoreS not found");
  auto *coreTubs = dynamic_cast<G4Tubs *>(coreVol->GetSolid());
  fR_Core     = coreTubs->GetOuterRadius();
  fFiberHalfZ = coreTubs->GetZHalfLength();

  std::cout << "Fiber optics (from material tables):\n"
            << "  S-fiber:       n_core=" << fN_CoreS << "  n_clad=" << fN_CladS
            << "  abs_len=" << fAbsLen_CoreS/m << " m\n"
            << "  C-fiber(Plastic): n_core=" << fN_CoreC << "  n_clad=" << fN_CladC
            << "  abs_len=" << fAbsLen_CoreC/m << " m\n"
            << "  C-fiber(Quartz):  n_core=" << fN_CoreQ << "  n_clad=" << fN_CladQ
            << "  abs_len=" << fAbsLen_CoreQ/m << " m\n"
            << "  R_core=" << fR_Core/mm << " mm"
            << "  half-z=" << fFiberHalfZ/cm << " cm" << std::endl;

  fOpticsInitialised = true;
}

void CaloXSteppingAction::fillOPInfo(const G4Step *step, bool verbose)
{
  if (!fOpticsInitialised) initOptics();

  // Optical constants read from material tables on first call
  const double N_CORE_S     = fN_CoreS;
  const double N_CLAD_S     = fN_CladS;
  const double N_CORE_C     = fN_CoreC;
  const double N_CLAD_C     = fN_CladC;
  const double N_CORE_Q     = fN_CoreQ;
  const double N_CLAD_Q     = fN_CladQ;
  const double R_CORE       = fR_Core;
  const double FIBER_HALF_Z = fFiberHalfZ;

  G4Track      *track        = step->GetTrack();
  G4StepPoint  *preStepPoint = step->GetPreStepPoint();
  G4StepPoint  *postStepPoint = step->GetPostStepPoint();
  const G4int   trackID      = track->GetTrackID();

  // Identify fiber sub-volume
  bool isCoreS = false, isCoreC = false, isCoreQ = false;
  bool isCladS = false, isCladC = false, isCladQ = false;
  const G4String detname =
      track->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName();
  if      (detname == "fiberCoreS")       isCoreS = true;
  else if (detname == "fiberCorePlastic") isCoreC = true;
  else if (detname == "fiberCoreQuartz")  isCoreQ = true;
  else if (detname == "fiberCladS")       isCladS = true;
  else if (detname == "fiberCladPlastic") isCladC = true;
  else if (detname == "fiberCladQuartz")  isCladQ = true;
  else { track->SetTrackStatus(fStopAndKill); return; }

  const int fiberIdx  = (isCladS || isCladC || isCladQ) ? 0 : 1;
  const int fiberNumber = preStepPoint->GetTouchableHandle()->GetCopyNumber(fiberIdx);
  const int rodNumber   = preStepPoint->GetTouchableHandle()->GetCopyNumber(fiberIdx + 2);
  const int layerNumber = preStepPoint->GetTouchableHandle()->GetCopyNumber(fiberIdx + 3);

  const bool isSampleRod = (rodNumber == opSampleRod && layerNumber == opSampleLayer);

  // -----------------------------------------------------------------------
  //  At photon birth: record info and compute analytical TIR for all rods.
  //  For non-sample rods, kill immediately after — no further tracking needed.
  // -----------------------------------------------------------------------
  if (track->GetCurrentStepNumber() == 1)
  {
    // Count Cerenkov photons born in core (all rods)
    if (isCoreS || isCoreC || isCoreQ)
    {
      auto *cp = track->GetCreatorProcess();
      if (cp && cp->GetProcessName() == "Cerenkov")
      {
        int fiberType = isCoreS ? 0 : (isCoreC ? 1 : 2); // 0=scint, 1=plastic, 2=quartz
        hh->accumulateOPsCer(fiberType, 1);
      }
    }

    // For non-sample rods, only core photons can be guided — skip clad-born photons
    if (!isSampleRod && !(isCoreS || isCoreC || isCoreQ))
    {
      track->SetTrackStatus(fStopAndKill);
      return;
    }

    // Build the photon record
    CaloXPhotonInfo photon;
    photon.trackID            = trackID;
    photon.productionPosition = preStepPoint->GetPosition() / cm;
    photon.productionMomentum = track->GetMomentum() / GeV;
    photon.productionTime     = track->GetGlobalTime() / ns;
    photon.polarization       = track->GetPolarization();
    photon.isCoreS = isCoreS; photon.isCoreC = isCoreC; photon.isCoreQ = isCoreQ;
    photon.isCladS = isCladS; photon.isCladC = isCladC; photon.isCladQ = isCladQ;
    photon.productionFiber    = fiberNumber;
    photon.productionRod      = rodNumber;
    photon.productionLayer    = layerNumber;
    auto *cp = track->GetCreatorProcess();
    if (cp)
    {
      photon.isCerenkov      = (cp->GetProcessName() == "Cerenkov");
      photon.isScintillation = (cp->GetProcessName() == "Scintillation");
    }

    // Analytical TIR for core photons (same logic, all rods)
    if (isCoreS || isCoreC || isCoreQ)
    {
      const double n_core = isCoreS ? N_CORE_S : (isCoreC ? N_CORE_C : N_CORE_Q);
      const double n_clad = isCoreS ? N_CLAD_S : (isCoreC ? N_CLAD_C : N_CLAD_Q);
      const double abs_len = isCoreS ? fAbsLen_CoreS : (isCoreC ? fAbsLen_CoreC : fAbsLen_CoreQ);
      computeAnalyticalTIR(photon,
                           track->GetMomentumDirection(),
                           preStepPoint->GetPosition(),
                           n_core, n_clad, R_CORE, FIBER_HALF_Z,
                           track->GetGlobalTime(), abs_len);
    }

    hh->photonData.push_back(photon);

    if (!isSampleRod)
    {
      track->SetTrackStatus(fStopAndKill);
      return;
    }
  }

  // -----------------------------------------------------------------------
  //  Sample rod only: update exit info when photon leaves the fiber
  // -----------------------------------------------------------------------
  if (postStepPoint->GetTouchableHandle()->GetVolume())
  {
    const G4String exitName =
        postStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
    if (exitName == "World" || exitName == "Calorimeter")
    {
      for (auto &photon : hh->photonData)
      {
        if (photon.trackID != trackID) continue;
        photon.exitPosition = postStepPoint->GetPosition() / cm;
        photon.exitMomentum = track->GetMomentum() / GeV;
        photon.exitTime     = track->GetGlobalTime() / ns;
        if (verbose)
          std::cout << "Sample-rod photon exit z="
                    << preStepPoint->GetPosition().z() / cm << " cm" << std::endl;
        break;
      }
    }
  }
}

