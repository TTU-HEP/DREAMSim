#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>

using namespace Pythia8;

// Simple method to do the filling of partons into the event record.

void fillPartons(int type, double ee, Event &event, ParticleData &pdt,
                 Rndm &rndm)
{

  // Reset event record to allow for new event.
  event.reset();

  // Information on a q qbar system, to be hadronized.
  if (type == 1)
  {
    int id = 2;
    double mm = pdt.m0(id);
    double pp = sqrtpos(ee * ee - mm * mm);
    event.append(id, 23, 101, 0, 0., 0., pp, ee, mm);
    event.append(-id, 23, 0, 101, 0., 0., -pp, ee, mm);

    // Information on a g g system, to be hadronized.
  }
  else if (type == 2)
  {
    event.append(21, 23, 101, 102, 0., 0., ee, ee);
    event.append(21, 23, 102, 101, 0., 0., -ee, ee);

    // Information on a g g g system, to be hadronized.
  }
}

//==========================================================================

int main()
{

  // Loop over kind of events to generate:
  // 1 = q qbar.
  // 2 = g g.
  int type = 1;

  // Set typical energy per parton.
  double ee = 50.0;

  // Set number of events to generate and to list.
  int nEvent = 10000;
  int nList = 1;

  // Generator; shorthand for event and particleData.
  Pythia pythia;
  Event &event = pythia.event;
  ParticleData &pdt = pythia.particleData;

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // Avoid the standard scrutiny of mother/daughter relations.
  // But note that other event checks are done below.
  pythia.readString("Check:event = off");

  // Optionally switch off resonance decays.
  // pythia.readString("ProcessLevel:resonanceDecays = off");

  // Optionally switch off ordinary decays.
  // pythia.readString("HadronLevel:Decay = off");

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Set true to also see space-time information in event listings.
  bool showScaleAndVertex = false;

  // Initialize.
  cout << "\n Now begin type = " << type << endl;
  // If Pythia fails to initialize, exit with error.
  if (!pythia.init())
    return 1;

  TFile *outFile = new TFile("PSParticles.root", "RECREATE");
  TTree *tree = new TTree("Particles", "final state particles");

  // output variables
  std::vector<int> particle_ID;
  std::vector<double> particle_E;
  std::vector<double> particle_px;
  std::vector<double> particle_py;
  std::vector<double> particle_pz;

  tree->Branch("particle_ID", &particle_ID);
  tree->Branch("particle_E", &particle_E);
  tree->Branch("particle_px", &particle_px);
  tree->Branch("particle_py", &particle_py);
  tree->Branch("particle_pz", &particle_pz);

  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {

    fillPartons(type, ee, event, pdt, pythia.rndm);

    // Generate events. Quit if failure.
    if (!pythia.next())
    {
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // List first few events.
    if (iEvent < nList)
    {
      event.list(showScaleAndVertex);
      // Also list junctions.
      event.listJunctions();
    }

    particle_ID.clear();
    particle_E.clear();
    particle_px.clear();
    particle_py.clear();
    particle_pz.clear();

    // Loop over all particles.
    for (int i = 0; i < event.size(); ++i)
    {
      const Particle &p = event[i];
      if (!p.isFinal())
        continue; // Only consider final state particles.
      if (p.pz() < 0.0)
        continue; // Only consider particles with positive pz.
      particle_ID.push_back(p.id());
      particle_E.push_back(p.e());
      particle_px.push_back(p.px());
      particle_py.push_back(p.py());
      particle_pz.push_back(p.pz());
    }

    tree->Fill(); // Fill the tree with the current event's data.

    // End of event loop.
  }

  tree->Write();    // Write the tree to the file.
  outFile->Close(); // Close the output file.
  delete outFile;   // Clean up the file pointer.

  // Print statistics and histograms.
  pythia.stat();

  // End of type loop and done.
  return 0;
}
