#include "globals.hh"
#include "G4ThreeVector.hh"

struct CaloXPhotonInfo
{
    G4int trackID = -1;
    G4ThreeVector productionPosition = G4ThreeVector(0, 0, 0);
    G4ThreeVector exitPosition = G4ThreeVector(0, 0, 0);
    G4ThreeVector productionMomentum = G4ThreeVector(0, 0, 0);
    G4ThreeVector exitMomentum = G4ThreeVector(0, 0, 0);
    G4ThreeVector polarization = G4ThreeVector(0, 0, 0);
    G4double productionTime = 0;
    G4double exitTime = 0;
    G4int productionFiber = -99;
    G4int productionRod   = -99;
    G4int productionLayer = -99;
    G4int exitFiber = -99;
    G4bool isCerenkov = false;
    G4bool isScintillation = false;
    G4bool isCoreC = false;
    G4bool isCoreS = false;
    G4bool isCladC = false;
    G4bool isCladS = false;

    // Analytical fast-tracking result (all fibers outside the sample rod)
    // Meridional approximation (pz only — exact only for rays through fiber axis)
    G4bool   isCaptured_m            = false; // TIR using meridional formula: |pz| >= n_clad/n_core
    G4bool   isAttenuated_m          = false; // true if killed by Beer-Lambert along meridional path
    G4double captureAngle_m          = 0.0;   // acos(|pz|) — angle with fiber axis (rad)
    G4double analyticalArrivalTime_m = 0.0;   // arrival time using pz only (ns)

    // Skew-ray (correct) calculation using conserved z-angular momentum L_z
    G4bool   isCaptured_s            = false; // TIR using sin(theta_wall) = sqrt(pz^2 + Lz^2/R^2)
    G4bool   isAttenuated_s          = false; // true if killed by Beer-Lambert along skew path
    G4double captureAngle_s          = 0.0;   // alpha = acos(|pz|), same definition as _m
    G4double analyticalArrivalTime_s = 0.0;   // arrival time using pz (ns)
};
