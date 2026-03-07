/// \file exampleCaloX.cc
/// \brief Main program of the CaloX example

// #include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
// #include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "QBBC.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
// #include "G4Cerenkov.hh"
#include "Randomize.hh"
#include <chrono>

#include "G4HadronicProcess.hh"
#include "G4GammaGeneralProcess.hh"

#include "CaloXDetectorConstruction.hh"
#include "CaloXPrimaryGeneratorAction.hh"
#include "CaloXRunAction.hh"
#include "CaloXEventAction.hh"
#include "CaloXSteppingAction.hh"
#include "CaloXTrackingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "CaloXTree.h"
#include "CaloXPhysicsList.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{

  bool batchJob = false;

  string macro;

  for (G4int i = 1; i < argc; i = i + 2)
  {
    string a = argv[i];
    if (G4String(argv[i]) == "-b")
    {
      macro = argv[i + 1];
      batchJob = true;
    }
    else if (G4String(argv[i]) == "-i")
    {
      macro = argv[i + 1];
      batchJob = false;
    }
    else if (a.substr(0, 1) != "-")
    {
      std::cout << "argument error: parameter shoudl start with -. " << a << std::endl;
      return 0;
    }
  }

  CaloXTree *histo = new CaloXTree(macro, argc, argv);

  G4UIExecutive *ui = nullptr;
  if (!batchJob)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Generate time-based random seeds so each run is statistically independent
  auto now = std::chrono::system_clock::now();
  auto micros = std::chrono::duration_cast<std::chrono::microseconds>(
                    now.time_since_epoch())
                    .count();
  long long t1 = micros / 1000000LL;
  long long t2 = micros % 1000000LL;
  std::cout << "Time-based seed components: t1=" << t1 << "   t2=" << t2 << std::endl;
  long seeds[2];
  int kseed = histo->getParamI("runNumber") + histo->getParamI("runSeq") * 3333;
  seeds[0] = long(t2) + long(kseed);
  seeds[1] = seeds[0] + 8134;

  // seeds[0]=2345;
  // seeds[1]=7999;

  G4Random::setTheSeeds(seeds);
  G4Random::showEngineStatus();
  std::cout << "seeds[0]=" << seeds[0] << "   seeds[1]=" << seeds[1] << std::endl;

  // Construct a serial run manager
  //
  // auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  auto *runManager = new G4RunManager;
  // auto *runManager = new G4MTRunManager;

  // Set mandatory initialization classes
  //
  auto detector = new CaloXDetectorConstruction(histo);
  runManager->SetUserInitialization(detector);

  //   optical physics from examples/extended/optical/OpNovice2
  //  auto physicsList = new FTFP_BERT;
  // Use QGSP_BERT as the base physics list.
  // CaloXPhysicsList (defined in include/CaloXPhysicsList.hh) extends QGSP_BERT
  // by removing the photonNuclear process; switch to it if that suppression is needed.
  auto physicsList = new QGSP_BERT;
  // auto physicsList = new CaloXPhysicsList();

  // physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics *opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);

  runManager->SetUserInitialization(physicsList);

  // G4Cerenkov* theCerenkovProcess=new G4Cerenkov("Cerenkov");
  // theCerenkovProcess->SetTrackSecondariesFirst(true);
  // nt MaxNumberPhotons=300;
  // theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumberPhotons);
  // physicsList->RegisterPhysics(theCerenkovProcess);

  auto gen_action = new CaloXPrimaryGeneratorAction(detector, histo);
  runManager->SetUserAction(gen_action);

  auto run_action = new CaloXRunAction(histo);
  runManager->SetUserAction(run_action);
  //
  auto event_action = new CaloXEventAction(detector, gen_action, histo);
  runManager->SetUserAction(event_action);
  //
  auto stepping_action = new CaloXSteppingAction(event_action, histo);
  runManager->SetUserAction(stepping_action);
  //
  auto tracking_action = new CaloXTrackingAction();
  runManager->SetUserAction(tracking_action);

  // runManager->SetNumberOfThreads(4);
  runManager->Initialize();
  // auto actionInitialization = new CaloXActionInitialization(detConstruction);
  // runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if (batchJob)
  {
    // batch mode
    std::cout << "debug:  batch mode" << std::endl;
    G4String command = "/control/execute ";
    cout << "command: " << command << endl;
    UImanager->ApplyCommand(command + macro); // macro does not have /run/beamOn 100

    //    beams are now defined in CaloXPrimaryGeneratorAction...
    // command="/gun/particle "+histo->getParamS("gun_particle");
    // cout<<"command: "<<command<<endl;
    // UImanager->ApplyCommand(command);

    // command="/gun/energy "+histo->getParamS("gun_energy")+" GeV";
    // cout<<"command: "<<command<<endl;
    // UImanager->ApplyCommand(command);;

    // string evtmax="100";
    command = "/run/beamOn " + histo->getParamS("numberOfEvents");
    cout << "command: " << command << endl;
    UImanager->ApplyCommand(command);
  }
  else
  {
    // interactive mode
    std::cout << "debug:  interactive mode" << std::endl;
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    std::cout << "debug:  interactive mode, step 2" << std::endl;
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  std::cout << "Job termination..." << std::endl;
  histo->EndJob();
  std::cout << "after histo->EndJob()..." << std::endl;

  G4Random::showEngineStatus();
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete histo;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
