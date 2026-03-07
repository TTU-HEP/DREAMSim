/// \file CaloXDetectorConstruction.cc
/// \brief Implementation of the CaloXDetectorConstruction class

#include "CaloXDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4MaterialPropertiesTable.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

#include "CaloXTree.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger *CaloXDetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloXDetectorConstruction::CaloXDetectorConstruction(CaloXTree *histo)
    : G4VUserDetectorConstruction(),
      hh(histo),
      fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CaloXDetectorConstruction::~CaloXDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *CaloXDetectorConstruction::Construct()
{
    // Define materials
    DefineMaterials();

    // Define volumes
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloXDetectorConstruction::DefineMaterials()
{
    // Lead material defined using NIST Manager
    auto nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_Fe");
    nistManager->FindOrBuildMaterial("G4_Cu");
    nistManager->FindOrBuildMaterial("G4_Pb");
    nistManager->FindOrBuildMaterial("G4_Si");
    nistManager->FindOrBuildMaterial("G4_W");
    nistManager->FindOrBuildMaterial("G4_PbWO4");
    nistManager->FindOrBuildMaterial("G4_BRASS");
    nistManager->FindOrBuildMaterial("G4_U");
    nistManager->FindOrBuildMaterial("G4_AIR");


    // (PolyVinylToluene, C_9H_10)
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    // H_Scintillator.
    auto mat_H = nistManager->FindOrBuildMaterial("G4_H");
    auto mat_C = nistManager->FindOrBuildMaterial("G4_C");
    double densitySC = 1.032 * g / cm3;
    G4Material *h_scintillator = new G4Material("H_Scintillator", densitySC, 2);
    h_scintillator->AddMaterial(mat_C, 0.91512109);
    h_scintillator->AddMaterial(mat_H, 0.084878906);

    // Liquid argon material
    G4double a; // mass of a mole;
    G4double z; // z=mean number of protons;
    G4double density;

    // Vacuum
    new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = universe_mean_density,
                   kStateGas, 2.73 * kelvin, 3.e-18 * pascal);

    // Print materials
    // std::cout << *(G4Material::GetMaterialTable()) << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *CaloXDetectorConstruction::DefineVolumes()
{
    // Geometry structure
    //   World                  Air
    //     - Calo               Air
    //        - Layer [80]      Cu
    //           - Rod [90]
    //               - Hole  [1]  Air
    //                  - CF [1]     Cherenkov Fiber
    //                  - SF [1]
    // data:
    //   cahnnel count:  80*100*(5+100+100) = 1.6 M
    //   timeT= (200.0-Zcalo)/20.0 + TOF
    //
    //   (event)
    //   nts  number of time slices
    //
    //   (hit)
    //   ID:  layerN*100+RodN
    //   edeprod  in rod (Cu)
    //   edepsc   in S-fiber
    //   edepch   in C-fiber
    //   sc
    //   scts[100]
    //   ch
    //   chts[100]
    //
    // Geometry parameters
    double fiberLength = 250.0 * cm;
    double holeDiameter = 0.25 * cm;
    double rodSize = 0.4 * cm;
    double noLayers = 80.0;
    double layerThickness = rodSize;
    double noRods = 90.0;

    double calorSizeX = rodSize * noRods;
    double calorSizeY = rodSize * noLayers;
    double calorSizeZ = fiberLength;

    // Use the calorimeter's space diagonal as the world half-size so that
    // any rotation (including 90° around Y) is always fully contained.
    double calorDiag  = std::sqrt(calorSizeX*calorSizeX +
                                  calorSizeY*calorSizeY +
                                  calorSizeZ*calorSizeZ);
    double worldSizeX = 1.2 * calorDiag;
    double worldSizeY = 1.2 * calorDiag;
    double worldSizeZ = 1.2 * calorDiag;

    //
    // World
    //
    G4Material *defaultMaterial = G4Material::GetMaterial("G4_AIR"); // G4_AIR or G4_Galactic

    G4VSolid *worldS = new G4Box("World",                                             // its name
                                 worldSizeX / 2.0, worldSizeY / 2.0, worldSizeZ / 2); // its size

    G4LogicalVolume *worldLV = new G4LogicalVolume(
        worldS,          // its solid
        defaultMaterial, // its material
        "World");        // its name

    G4VPhysicalVolume *worldPV = new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        worldLV,         // its logical volume
        "World",         // its name
        0,               // its mother volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    //
    // Calorimeter
    //
    auto calorMaterial = G4Material::GetMaterial("G4_Cu"); // CaloX nominal
                                                           // auto calorMaterial = G4Material::GetMaterial("G4_PbWO4");
                                                           // auto calorMaterial = G4Material::GetMaterial("G4_Si");
                                                           // auto calorMaterial = G4Material::GetMaterial("G4_W");
                                                           // auto calorMaterial = G4Material::GetMaterial("G4_Pb");
                                                           // auto calorMaterial = G4Material::GetMaterial("G4_U");
                                                           // auto sensorMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
                                                           // auto sensorMaterial = G4Material::GetMaterial("H_Scintillator");

    G4NistManager *nistManager = G4NistManager::Instance();

    G4int ncomponents, natoms;
    // Elements for fiber materials
    G4Element *H = nistManager->FindOrBuildElement(1);
    G4Element *C = nistManager->FindOrBuildElement(6);
    G4Element *N = nistManager->FindOrBuildElement(7);
    G4Element *O = nistManager->FindOrBuildElement(8);
    G4Element *F = nistManager->FindOrBuildElement(9);
    G4Element *Si = nistManager->FindOrBuildElement(14);

    auto calorimeterS = new G4Box("Calorimeter",                                      // its name
                                  calorSizeX / 2., calorSizeY / 2., calorSizeZ / 2.); // its size

    auto calorLV = new G4LogicalVolume(
        calorimeterS,   // its solid
        calorMaterial,  // its material
        "Calorimeter"); // its name

    G4RotationMatrix *xRot = new G4RotationMatrix; // Rotates X and Z axes only
    xRot->rotateX(hh->getParamF("caloRotationX") * deg);
    xRot->rotateY(hh->getParamF("caloRotationY") * deg);
    xRot->rotateZ(0. * deg);

    new G4PVPlacement(
        xRot,            // rotate by caloRotationX/Y (2') degree
        G4ThreeVector(), // at (0,0,0)
        calorLV,         // its logical volume
        "Calorimeter",   // its name
        worldLV,         // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    //
    // Layer
    //
    auto layerS = new G4Box("Layer",                                                  // its name
                            calorSizeX / 2.0, layerThickness / 2.0, calorSizeZ / 2.); // its size

    auto layerLV = new G4LogicalVolume(
        layerS,        // its solid
        calorMaterial, // its material
        "Layer");      // its name

    new G4PVReplica(
        "Layer",         // its name
        layerLV,         // its logical volume
        calorLV,         // its mother
        kYAxis,          // axis of replication
        noLayers,        // number of replic
        layerThickness); // width of replica

    //
    // Rods in a layer
    //
    auto rodS = new G4Box("Rod",                                          // its name
                          rodSize / 2.0, rodSize / 2.0, calorSizeZ / 2.); // its size

    auto rodLV = new G4LogicalVolume(
        rodS,          // its solid  [FIXED: was layerS — wrong geometry]
        calorMaterial, // its material
        "Rod");        // its name

    new G4PVReplica(
        "Rod",    // its name
        rodLV,    // its logical volume
        layerLV,  // its mother
        kXAxis,   // axis of replication
        noRods,   // number of replic
        rodSize); // witdth of replica

    //
    // Hole in a Rod
    //

    G4Material *holeMaterial = G4Material::GetMaterial("G4_AIR"); // G4_AIR or G4_Galactic

    auto holeS = new G4Tubs("Hole", // its name
                            0.0, holeDiameter / 2.0, calorSizeZ / 2., 0.0 * deg, 360. * deg);

    auto holeLV = new G4LogicalVolume(
        holeS,        // its solid
        holeMaterial, // its material
        "Hole");      // its name

    new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        holeLV,          // its logical volume
        "Hole",          // its name
        rodLV,           // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    //
    //  Fibers
    //

    ///--- for scintillation fiber core ---
    double density;
    G4Material *polystyrene =
        new G4Material("Polystyrene", density = 1.05 * g / cm3, ncomponents = 2);
    polystyrene->AddElement(C, natoms = 8);
    polystyrene->AddElement(H, natoms = 8);

    ///--- for cladding (scintillation fibers) ---
    G4Material *pmma_clad =
        new G4Material("PMMA_Clad", density = 1.19 * g / cm3, ncomponents = 3);
    pmma_clad->AddElement(C, natoms = 5);
    pmma_clad->AddElement(H, natoms = 8);
    pmma_clad->AddElement(O, natoms = 2);

    ///--- for Cherenkov fiber core ---
    G4Material *pmma =
        new G4Material("PMMA", density = 1.19 * g / cm3, ncomponents = 3);
    pmma->AddElement(C, natoms = 5);
    pmma->AddElement(H, natoms = 8);
    pmma->AddElement(O, natoms = 2);

    ///--- for cladding (Cerenkov fibers) ---
    G4Material *fluorinatedPolymer =
        new G4Material("Fluorinated_Polymer", density = 1.43 * g / cm3, ncomponents = 2);
    fluorinatedPolymer->AddElement(C, 2);
    fluorinatedPolymer->AddElement(F, 2);
    G4MaterialPropertiesTable *mpPMMA;
    G4MaterialPropertiesTable *mpFS;
    G4MaterialPropertiesTable *mpPS;

    //--- Generate and add material properties table ---
    G4double PhotonEnergy[] = {2.00 * eV, 2.03 * eV, 2.06 * eV, 2.09 * eV, 2.12 * eV,
                               2.15 * eV, 2.18 * eV, 2.21 * eV, 2.24 * eV, 2.27 * eV,
                               2.30 * eV, 2.33 * eV, 2.36 * eV, 2.39 * eV, 2.42 * eV,
                               2.45 * eV, 2.48 * eV, 2.51 * eV, 2.54 * eV, 2.57 * eV,
                               2.60 * eV, 2.63 * eV, 2.66 * eV, 2.69 * eV, 2.72 * eV,
                               2.75 * eV, 2.78 * eV, 2.81 * eV, 2.84 * eV, 2.87 * eV,
                               2.90 * eV, 2.93 * eV, 2.96 * eV, 2.99 * eV, 3.02 * eV,
                               3.05 * eV, 3.08 * eV, 3.11 * eV, 3.14 * eV, 3.17 * eV,
                               3.20 * eV, 3.23 * eV, 3.26 * eV, 3.29 * eV, 3.32 * eV,
                               3.35 * eV, 3.38 * eV, 3.41 * eV, 3.44 * eV, 3.47 * eV};

    const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

    //--- PMMA ---
    G4double RefractiveIndex_PMMA[nEntries] =
        {
            1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468,
            1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468,
            1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468,
            1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468,
            1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468, 1.468};
    mpPMMA = new G4MaterialPropertiesTable();
    mpPMMA->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_PMMA, nEntries);

    G4double Absorption_PMMA[nEntries];
    std::fill_n(Absorption_PMMA, nEntries, 5.0 * m);
    mpPMMA->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_PMMA, nEntries);

    pmma->SetMaterialPropertiesTable(mpPMMA);
    pmma_clad->SetMaterialPropertiesTable(mpPMMA);

    //--- Fluorinated Polymer (FS) ---
    G4double RefractiveIndex_FluorinatedPolymer[nEntries] =
        {
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42};
    mpFS = new G4MaterialPropertiesTable();
    mpFS->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_FluorinatedPolymer, nEntries);

    G4double Absorption_FluorinatedPolymer[nEntries];
    std::fill_n(Absorption_FluorinatedPolymer, nEntries, 10.0 * m);
    mpFS->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_FluorinatedPolymer, nEntries);

    fluorinatedPolymer->SetMaterialPropertiesTable(mpFS);

    // -- Polystyrene --
    G4double RefractiveIndex_Polystyrene[nEntries] =
        {
            1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
            1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
            1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
            1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
            1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59};
    mpPS = new G4MaterialPropertiesTable();
    mpPS->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Polystyrene, nEntries);
    G4double Absorption_Polystyrene[nEntries];
    std::fill_n(Absorption_Polystyrene, nEntries, 3.0 * m);
    mpPS->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_Polystyrene, nEntries);

    polystyrene->SetMaterialPropertiesTable(mpPS);

    //---Materials for Cerenkov fiber---
    G4Material *clad_C_Material = fluorinatedPolymer;
    G4Material *core_C_Material = pmma;
    //---Materials for Scintillation fiber---
    G4Material *clad_S_Material = pmma_clad;
    G4Material *core_S_Material = polystyrene;

    // Parameters for fibers
    double clad_C_rMin = 0.39 * mm; // cladding cherenkov minimum radius
    double clad_C_rMax = 0.40 * mm; // cladding cherenkov max radius
    // double clad_C_Dz = fiberLength / 2.0; // cladding cherenkov lenght
    // double clad_C_Sphi = 0.;              // cladding cherenkov min rotation
    // double clad_C_Dphi = 2. * M_PI;       // cladding chrenkov max rotation

    double core_C_rMin = 0. * mm;
    double core_C_rMax = 0.39 * mm;
    // double core_C_Dz = clad_C_Dz;
    // double core_C_Sphi = 0.;
    // double core_C_Dphi = 2. * M_PI;

    double clad_S_rMin = 0.39 * mm;
    double clad_S_rMax = 0.40 * mm;
    // double clad_S_Dz = clad_C_Dz;
    // double clad_S_Sphi = 0.;
    // double clad_S_Dphi = 2. * M_PI;

    double core_S_rMin = 0. * mm;
    double core_S_rMax = 0.39 * mm;
    // double core_S_Dz = clad_C_Dz;
    // double core_S_Sphi = 0.;
    // double core_S_Dphi = 2. * M_PI;

    // double theta_unit = 0;
    // double deltatheta = 0;
    // double thetaofcenter = 0;

    // creating fibers solids
    // G4cout << "r_clad= " << clad_C_rMax << " r_coreC=" << core_C_rMax << " r_coreS=" << core_S_rMax << G4endl;
    auto fiber = new G4Tubs("Fiber", 0, clad_C_rMax, fiberLength / 2., 0 * deg, 360. * deg); // S is the same
    auto fiberC = new G4Tubs("fiberC", 0, core_C_rMax, fiberLength / 2., 0 * deg, 360. * deg);
    auto fiberS = new G4Tubs("fiberS", 0, core_S_rMax, fiberLength / 2., 0 * deg, 360. * deg);

    auto fiberCLog = new G4LogicalVolume(fiber, clad_C_Material, "fiberCladC");
    auto fiberSLog = new G4LogicalVolume(fiber, clad_S_Material, "fiberCladS");

    G4LogicalVolume *fiberCoreCLog = new G4LogicalVolume(fiberC, core_C_Material, "fiberCoreC");
    G4LogicalVolume *fiberCoreSLog = new G4LogicalVolume(fiberS, core_S_Material, "fiberCoreS");

    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fiberCoreCLog, "fiberCoreCherePhys", fiberCLog, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fiberCoreSLog, "fiberCoreScintPhys", fiberSLog, false, 0);

    double R = clad_C_rMax * 2.0 + 0.01; // 10 micron gap between cenral and peripheral fibers
    double cx1 = R * cos(30.0 * deg);
    double cy1 = R * sin(30.0 * deg);
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), fiberCLog, "fiberCladC", holeLV, false, 0, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(cx1, cy1, 0.0), fiberCLog, "fiberCladC", holeLV, false, 1, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(-cx1, cy1, 0.0), fiberCLog, "fiberCladC", holeLV, false, 2, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0., -R, 0.), fiberCLog, "fiberCladC", holeLV, false, 3, fCheckOverlaps);

    new G4PVPlacement(0, G4ThreeVector(cx1, -cy1, 0.), fiberSLog, "fiberCladS", holeLV, false, 1, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0.0, R, 0.), fiberSLog, "fiberCladS", holeLV, false, 2, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(-cx1, -cy1, 0.), fiberSLog, "fiberCladS", holeLV, false, 3, fCheckOverlaps);

    /*if(sd){
     fiberCoreCLog->SetSensitiveDetector(sd);
     fiberCoreSLog->SetSensitiveDetector(sd);
     }*/

    //
    // print parameters
    //
    std::cout << "'' " << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "---> The calorimeter is " << noLayers << " layers of: [ " << std::endl;
    std::cout << layerThickness / mm << "mm of (layer) " << calorMaterial->GetName() << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    //
    // Visualization attributes
    //
    // worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

    worldLV->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.0, 0.0, 1.0, 0.5)));  // blue
    calorLV->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(1.0, 0.0, 0.0, 0.1)));  // red
    layerLV->SetVisAttributes(new G4VisAttributes(FALSE, G4Colour(0.0, 1.0, 0.0, 0.6))); // green
    rodLV->SetVisAttributes(new G4VisAttributes(FALSE, G4Colour(0.0, 0.0, 0.0, 0.6)));   // blue
    // holeLV->SetVisAttributes(new G4VisAttributes(FALSE,G4Colour(1.0,1.0,1.0))); // black
    holeLV->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(1.0, 1.0, 1.0, 0.5))); // white
    fiberCLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.8, 0.5, 0.8, 0.9)));
    fiberCoreCLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.98, 0.5, 0.98, 0.9)));
    fiberSLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.0, 0.5, 0.8, 0.9)));       // red
    fiberCoreSLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.0, 0.98, 0.98, 0.9))); // red

    //
    // Always return the physical World
    //
    return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CaloXDetectorConstruction::ConstructSDandField()
{
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue;
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
