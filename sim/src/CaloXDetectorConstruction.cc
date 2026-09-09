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
#include "CaloXFiberMap.h"
#include <map>
#include <sstream>
#include <cstdlib>

G4ThreadLocal G4GlobalMagFieldMessenger *CaloXDetectorConstruction::fMagFieldMessenger = nullptr;

CaloXDetectorConstruction::CaloXDetectorConstruction(CaloXTree *histo)
    : G4VUserDetectorConstruction(),
      hh(histo),
      fCheckOverlaps(true)
{
    //  The fiber map says what each copper contains.  It is required: building
    //  the calorimeter without it would silently give a detector that is not the
    //  one being modelled, so stop instead of guessing.
    const std::string mapFile = hh->getParamS("fiberMapFile", false, "data/fibermap.json");
    if (!fFiberMap.load(mapFile))
    {
        std::cout << "CaloXDetectorConstruction: could not read the fiber map from '"
                  << mapFile << "'.  Set fiberMapFile in the mac file (or pass "
                  << "-fiberMapFile <path>) to point at it.  Exiting." << std::endl;
        std::exit(1);
    }
    std::cout << "CaloXDetectorConstruction: fiber map from " << mapFile << std::endl;
    fFiberMap.print();
}

CaloXDetectorConstruction::~CaloXDetectorConstruction()
{
}

G4VPhysicalVolume *CaloXDetectorConstruction::Construct()
{
    // Define materials
    DefineMaterials();

    // Define volumes
    return DefineVolumes();
}

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

    // std::cout << *(G4Material::GetMaterialTable()) << std::endl;
}

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
    //  Taken from the fiber map so the two can never disagree (nominally 80 x 90).
    double noLayers = fFiberMap.nLayers();
    double layerThickness = rodSize;
    double noRods = fFiberMap.nRods();

    double calorSizeX = rodSize * noRods;
    double calorSizeY = rodSize * noLayers;
    double calorSizeZ = fiberLength;

    // Use the calorimeter's space diagonal as the world half-size so that
    // any rotation (including 90° around Y) is always fully contained.
    double calorDiag = std::sqrt(calorSizeX * calorSizeX +
                                 calorSizeY * calorSizeY +
                                 calorSizeZ * calorSizeZ);
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

    //  Layers, rods and holes are built after the fibers, further down: what a
    //  copper contains depends on the fiber map, so the fiber logical volumes
    //  have to exist first.

    //
    //  Fibers
    //
    // Fiber material summary:
    //   Fiber      Component   Material              Elements   Density (g/cm3)   n       Att. length
    //   ---------  ---------   --------------------  ---------  ---------------   ------  -----------
    //   S-fiber    core        Polystyrene           C8H8       1.05              1.622   2.0 m
    //              clad        PMMA_Clad             C5H8O2     1.19              1.504   5.0 m
    //   C-Plastic  core        PMMA                  C5H8O2     1.19              1.504   5.0 m
    //              clad        Fluorinated_Polymer   C2F2       1.43              1.42    10.0 m
    //   C-Quartz   core        Fused_Silica          SiO2       2.2              1.468   10.0 m
    //              clad        Hard_Polymer          C2F2       1.43              1.42    10.0 m

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

    ///--- for Cherenkov Quartz fiber core (Fused Silica, SiO2) ---
    G4Material *fusedSilica =
        new G4Material("Fused_Silica", density = 2.2 * g / cm3, ncomponents = 2);
    fusedSilica->AddElement(Si, natoms = 1);
    fusedSilica->AddElement(O, natoms = 2);

    ///--- for Cherenkov Quartz fiber cladding (Hard Polymer) ---
    G4Material *hardPolymer =
        new G4Material("Hard_Polymer", density = 1.43 * g / cm3, ncomponents = 2);
    hardPolymer->AddElement(C, 2);
    hardPolymer->AddElement(F, 2);

    G4MaterialPropertiesTable *mpPMMA;
    G4MaterialPropertiesTable *mpFS;
    G4MaterialPropertiesTable *mpPS;
    G4MaterialPropertiesTable *mpFusedSilica;
    G4MaterialPropertiesTable *mpHardPolymer;

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
            1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504,
            1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504,
            1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504,
            1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504,
            1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504, 1.504};
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
            1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622,
            1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622,
            1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622,
            1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622,
            1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622, 1.622};
    mpPS = new G4MaterialPropertiesTable();
    mpPS->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Polystyrene, nEntries);
    G4double Absorption_Polystyrene[nEntries];
    std::fill_n(Absorption_Polystyrene, nEntries, 2.0 * m);
    mpPS->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_Polystyrene, nEntries);

    polystyrene->SetMaterialPropertiesTable(mpPS);

    //--- Fused Silica (Cherenkov Quartz fiber core) ---
    G4double RefractiveIndex_FusedSilica[nEntries];
    std::fill_n(RefractiveIndex_FusedSilica, nEntries, 1.468);
    mpFusedSilica = new G4MaterialPropertiesTable();
    mpFusedSilica->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_FusedSilica, nEntries);
    G4double Absorption_FusedSilica[nEntries];
    std::fill_n(Absorption_FusedSilica, nEntries, 10.0 * m);
    mpFusedSilica->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_FusedSilica, nEntries);
    fusedSilica->SetMaterialPropertiesTable(mpFusedSilica);

    //--- Hard Polymer (Cherenkov Quartz fiber cladding) ---
    G4double RefractiveIndex_HardPolymer[nEntries];
    std::fill_n(RefractiveIndex_HardPolymer, nEntries, 1.42);
    mpHardPolymer = new G4MaterialPropertiesTable();
    mpHardPolymer->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_HardPolymer, nEntries);
    G4double Absorption_HardPolymer[nEntries];
    std::fill_n(Absorption_HardPolymer, nEntries, 10.0 * m);
    mpHardPolymer->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_HardPolymer, nEntries);
    hardPolymer->SetMaterialPropertiesTable(mpHardPolymer);

    //---Materials for Cerenkov Plastic fiber---
    G4Material *clad_Plastic_Material = fluorinatedPolymer;
    G4Material *core_Plastic_Material = pmma;
    //---Materials for Cerenkov Quartz fiber (Fused Silica core, Hard Polymer cladding)---
    G4Material *clad_Quartz_Material = hardPolymer;
    G4Material *core_Quartz_Material = fusedSilica;
    //---Materials for Scintillation fiber---
    G4Material *clad_S_Material = pmma_clad;
    G4Material *core_S_Material = polystyrene;

    // Parameters for fibers
    double clad_Plastic_rMin = 0.39 * mm; // cladding cherenkov minimum radius
    double clad_Plastic_rMax = 0.40 * mm; // cladding cherenkov max radius
    // double clad_C_Dz = fiberLength / 2.0; // cladding cherenkov lenght
    // double clad_C_Sphi = 0.;              // cladding cherenkov min rotation
    // double clad_C_Dphi = 2. * M_PI;       // cladding chrenkov max rotation

    double core_Plastic_rMin = 0. * mm;
    double core_Plastic_rMax = 0.39 * mm;
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
    // G4cout << "r_clad= " << clad_Plastic_rMax << " r_coreC=" << core_Plastic_rMax << " r_coreS=" << core_S_rMax << G4endl;
    auto fiber = new G4Tubs("Fiber", 0, clad_Plastic_rMax, fiberLength / 2., 0 * deg, 360. * deg); // S is the same
    auto fiberC = new G4Tubs("fiberC", 0, core_Plastic_rMax, fiberLength / 2., 0 * deg, 360. * deg);
    auto fiberS = new G4Tubs("fiberS", 0, core_S_rMax, fiberLength / 2., 0 * deg, 360. * deg);

    auto fiberPlasticLog = new G4LogicalVolume(fiber, clad_Plastic_Material, "fiberCladPlastic");
    auto fiberQuartzLog = new G4LogicalVolume(fiber, clad_Quartz_Material, "fiberCladQuartz");
    auto fiberSLog = new G4LogicalVolume(fiber, clad_S_Material, "fiberCladS");

    G4LogicalVolume *fiberCorePlasticLog = new G4LogicalVolume(fiberC, core_Plastic_Material, "fiberCorePlastic");
    G4LogicalVolume *fiberCoreQuartzLog = new G4LogicalVolume(fiberC, core_Quartz_Material, "fiberCoreQuartz");
    G4LogicalVolume *fiberCoreSLog = new G4LogicalVolume(fiberS, core_S_Material, "fiberCoreS");

    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fiberCorePlasticLog, "fiberCoreCherePlasticPhys", fiberPlasticLog, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fiberCoreQuartzLog, "fiberCoreChereQuartzPhys", fiberQuartzLog, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fiberCoreSLog, "fiberCoreScintPhys", fiberSLog, false, 0);

    //  Seven fiber slots per copper: one on the axis and six at 30, 90, 150, 210,
    //  270 and 330 degrees.  Four of them carry Cherenkov fibers and three carry
    //  scintillating ones, and whether the Cherenkov four are quartz or plastic
    //  is what the fiber map decides per copper.
    //
    //  (The 150-degree slot used to be left empty while the third S fiber was
    //  placed at 210 degrees, exactly on top of a quartz fiber; navigation
    //  resolved to whichever was placed first, so that S fiber was invisible and
    //  each copper had 4 C against 2 effective S.  fCheckOverlaps cannot see
    //  this: G4 samples points on the surface of the new solid, which for two
    //  exactly coincident solids return kSurface, not kInside.)
    double R = clad_Plastic_rMax * 2.0 + 0.01; // 10 micron gap between cenral and peripheral fibers
    double cx1 = R * cos(30.0 * deg);
    double cy1 = R * sin(30.0 * deg);

    const G4ThreeVector cherenkovSlot[4] = {
        G4ThreeVector(0., 0., 0.),        //  on the axis
        G4ThreeVector(cx1, cy1, 0.),      //   30 degrees
        G4ThreeVector(-cx1, -cy1, 0.),    //  210 degrees
        G4ThreeVector(0., -R, 0.)};       //  270 degrees
    const G4ThreeVector scintSlot[3] = {
        G4ThreeVector(cx1, -cy1, 0.),     //  330 degrees
        G4ThreeVector(0., R, 0.),         //   90 degrees
        G4ThreeVector(-cx1, cy1, 0.)};    //  150 degrees

    G4Material *holeMaterial = G4Material::GetMaterial("G4_AIR"); // G4_AIR or G4_Galactic
    auto holeS = new G4Tubs("Hole", 0.0, holeDiameter / 2.0, calorSizeZ / 2.,
                            0.0 * deg, 360. * deg);

    //  One hole logical volume per Cherenkov flavour.  Copy numbers are kept as
    //  they were: Cherenkov fibers are 0..3 and scintillating fibers are 1..3.
    G4LogicalVolume *holeLV[2] = {nullptr, nullptr}; //  [0] plastic, [1] quartz
    for (int flavour = 0; flavour < 2; ++flavour)
    {
        const bool quartz = (flavour == 1);
        holeLV[flavour] = new G4LogicalVolume(
            holeS, holeMaterial, quartz ? "HoleQuartz" : "HolePlastic");
        G4LogicalVolume *cherenkovLV = quartz ? fiberQuartzLog : fiberPlasticLog;
        const G4String cherenkovName = quartz ? "fiberCladQuartz" : "fiberCladPlastic";
        for (int i = 0; i < 4; ++i)
            new G4PVPlacement(0, cherenkovSlot[i], cherenkovLV, cherenkovName,
                              holeLV[flavour], false, i, fCheckOverlaps);
        for (int i = 0; i < 3; ++i)
            new G4PVPlacement(0, scintSlot[i], fiberSLog, "fiberCladS",
                              holeLV[flavour], false, i + 1, fCheckOverlaps);
    }

    //
    //  Rods, and the layers that hold them
    //
    //  A rod is copper with at most one hole in it, so there are only three kinds:
    //  a hole with quartz Cherenkov fibers, a hole with plastic ones, and solid
    //  copper.  Every rod is placed explicitly rather than replicated, because
    //  replicas share one logical volume and so cannot differ from each other.
    //  The physical volume is called "Rod" in all three cases: the stepping
    //  action recognises rods by that name.
    //  Rods tile the calorimeter exactly (calorSizeX = noRods * rodSize, and the
    //  same in y), so the copper of the Calorimeter and Layer volumes is never
    //  actually reached: what a rod is made of is what is there.  A rod outside
    //  the detector outline is therefore filled with air rather than copper, and
    //  is called "AirGap" instead of "Rod" so that the stepping action does not
    //  count it as absorber -- it falls through to caloType 0, which contributes
    //  to eCalotruth but not to eRodtruth.
    auto rodS = new G4Box("Rod", rodSize / 2.0, rodSize / 2.0, calorSizeZ / 2.);
    G4LogicalVolume *rodLV[CaloXFiberMap::kNTypes] = {nullptr, nullptr, nullptr, nullptr};
    const char *rodLogName[CaloXFiberMap::kNTypes] = {"RodEmpty", "RodPlastic", "RodQuartz",
                                                      "AirGap"};
    const char *rodPhysName[CaloXFiberMap::kNTypes] = {"Rod", "Rod", "Rod", "AirGap"};
    for (int t = 0; t < CaloXFiberMap::kNTypes; ++t)
    {
        const bool absent = (t == CaloXFiberMap::kAbsent);
        rodLV[t] = new G4LogicalVolume(rodS, absent ? defaultMaterial : calorMaterial,
                                       rodLogName[t]);
        if (t == CaloXFiberMap::kPlastic || t == CaloXFiberMap::kQuartz)
            new G4PVPlacement(0, G4ThreeVector(), holeLV[t == CaloXFiberMap::kQuartz ? 1 : 0],
                              "Hole", rodLV[t], false, 0, fCheckOverlaps);
    }

    auto layerS = new G4Box("Layer", calorSizeX / 2.0, layerThickness / 2.0, calorSizeZ / 2.);

    //  Layers that hold the same sequence of rod types can share one logical
    //  volume, which keeps the number of placements down: outside the central
    //  region four consecutive layers are identical by construction.
    std::map<std::string, G4LogicalVolume *> layerByPattern;
    const int nRodsI = int(noRods);
    const int nLayersI = int(noLayers);

    for (int layer = 0; layer < nLayersI; ++layer)
    {
        std::string pattern(nRodsI, '0');
        for (int rod = 0; rod < nRodsI; ++rod)
            pattern[rod] = char('0' + int(fFiberMap.type(rod, layer)));

        G4LogicalVolume *thisLayerLV = nullptr;
        std::map<std::string, G4LogicalVolume *>::iterator known = layerByPattern.find(pattern);
        if (known != layerByPattern.end())
        {
            thisLayerLV = known->second;
        }
        else
        {
            std::ostringstream lname;
            lname << "Layer" << layerByPattern.size();
            thisLayerLV = new G4LogicalVolume(layerS, calorMaterial, lname.str());
            for (int rod = 0; rod < nRodsI; ++rod)
            {
                const int t = pattern[rod] - '0';
                new G4PVPlacement(
                    0,
                    G4ThreeVector((rod + 0.5) * rodSize - calorSizeX / 2.0, 0., 0.),
                    rodLV[t], rodPhysName[t], thisLayerLV, false, rod, fCheckOverlaps);
            }
            layerByPattern[pattern] = thisLayerLV;
            thisLayerLV->SetVisAttributes(new G4VisAttributes(FALSE, G4Colour(0.0, 1.0, 0.0, 0.6)));
        }

        new G4PVPlacement(
            0,
            G4ThreeVector(0., (layer + 0.5) * layerThickness - calorSizeY / 2.0, 0.),
            thisLayerLV, "Layer", calorLV, false, layer, fCheckOverlaps);
    }

    std::cout << "CaloXDetectorConstruction: built " << nLayersI << " layers from "
              << layerByPattern.size() << " distinct rod patterns, "
              << layerByPattern.size() * nRodsI << " rod placements" << std::endl;

    /*if(sd){
     fiberCorePlasticLog->SetSensitiveDetector(sd);
     fiberCoreSLog->SetSensitiveDetector(sd);
     }*/

    //
    // Visualization attributes
    //
    // worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

    worldLV->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.0, 0.0, 1.0, 0.5)));  // blue
    calorLV->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(1.0, 0.0, 0.0, 0.1)));  // red
    for (int t = 0; t < CaloXFiberMap::kNTypes; ++t)
        rodLV[t]->SetVisAttributes(new G4VisAttributes(FALSE, G4Colour(0.0, 0.0, 0.0, 0.6))); // blue
    for (int flavour = 0; flavour < 2; ++flavour)
        holeLV[flavour]->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(1.0, 1.0, 1.0, 0.5))); // white
    fiberPlasticLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.8, 0.5, 0.8, 0.9)));
    fiberCorePlasticLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.98, 0.5, 0.98, 0.9)));
    fiberQuartzLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.5, 0.8, 0.5, 0.9)));
    fiberCoreQuartzLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.5, 0.98, 0.5, 0.9)));
    fiberSLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.0, 0.5, 0.8, 0.9)));       // red
    fiberCoreSLog->SetVisAttributes(new G4VisAttributes(TRUE, G4Colour(0.0, 0.98, 0.98, 0.9))); // red

    //
    // Always return the physical World
    //
    return worldPV;
}

void CaloXDetectorConstruction::ConstructSDandField()
{
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue;
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(0);

    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}
