/// \file CaloXSteppingAction.hh
/// \brief Definition of the CaloXSteppingAction class

#ifndef CaloXSteppingAction_h
#define CaloXSteppingAction_h 1

#include "G4UserSteppingAction.hh"

#include "CaloXEventAction.hh"
// class CaloXEventAction;

class CaloXID;
class CaloXHit;
class CaloXTree;

/// Stepping action class.
///
/// In UserSteppingAction() there are collected the energy deposit and track
/// lengths of charged particles in Absober and Gap layers and
/// updated in CaloXRunData object.

class CaloXSteppingAction : public G4UserSteppingAction
{
public:
  CaloXSteppingAction(CaloXEventAction *eventAction, CaloXTree *histo);
  virtual ~CaloXSteppingAction();

  virtual void UserSteppingAction(const G4Step *step);

private:
  CaloXEventAction *fEventAction;
  CaloXTree *hh;

  double getBirk(const G4Step *step);
  double getBirkHC(double dEStep, double step, double charge, double density);
  double getBirkL3(double dEStep, double step, double charge, double density);
  vector<double> UserCerenkov(const G4Step *step);

  //  SiPM PDE handling...
  void initPDE(int sipmType); // 1= J 6mm 6.0V,  2=J 6 mm 2.5V
  float getPDE(float lambda);

  int pdeLambdaMin, pdeLambdaMax;
  std::vector<float> sipmPDE;

  // Optical photon sampling rod/layer (configurable, defaults: rod=45, layer=40)
  int opSampleRod;
  int opSampleLayer;

  // Refractive indices, absorption lengths, and geometry —
  // populated lazily on first step (after Construct())
  bool   fOpticsInitialised = false;
  double fN_CoreS, fN_CladS;   // S-fiber: Polystyrene core, PMMA_Clad cladding
  double fN_CoreC, fN_CladC;   // Plastic C-fiber: PMMA core, Fluorinated_Polymer cladding
  double fN_CoreQ, fN_CladQ;   // Quartz C-fiber: Fused_Silica core, Hard_Polymer cladding
  double fAbsLen_CoreS;         // bulk absorption length, Polystyrene (mm)
  double fAbsLen_CoreC;         // bulk absorption length, PMMA (mm)
  double fAbsLen_CoreQ;         // bulk absorption length, Fused_Silica (mm)
  double fR_Core;               // core radius, read from fiberCoreS logical volume
  double fFiberHalfZ;           // fiber half-length, read from fiberCoreS logical volume
  void   initOptics();          // performs the one-time lookup
  double findInvisible(const G4Step *step, bool verbose = false);
  void fillOPInfo(const G4Step *step, bool verbose = false);
};


#endif
