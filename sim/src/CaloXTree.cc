#include "CaloXTree.h"

#include <chrono>  // from std::
#include <cstdlib> // for rand() on archer.
#include <ctime>
#include <fstream>  // for input/output files
#include <iomanip>  // for setw() in cout,
#include <iostream> // for cout
#include <math.h>   // for sin(x) etc.
#include <sstream>  // for string stream
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include <numeric>

#include "CaloXHit.h"
#include "CaloXID.h"
#include "CaloXPhotonInfo.h"

using namespace std;

// ------------------------------------------------------------------

CaloXTree::CaloXTree(string macFileName, int argc, char **argv)
{
  cout << "initializing CaloXTree...   macFileName:" << macFileName << endl;

  readMacFile(macFileName);

  //  overwrite params from argc, argv...
  for (int i = 1; i < argc; i = i + 2)
  {
    if (string(argv[i]) == "-b" || string(argv[i]) == "-i")
      continue;
    if (i + 1 >= argc) {
      std::cout << "CaloXTree: WARNING: flag '" << argv[i]
                << "' has no value — ignoring." << std::endl;
      break;
    }
    string a = argv[i];
    string b = argv[i + 1];
    setParam(a.substr(1, a.size() - 1), b);
  }

  // Print final effective parameter values after all overrides
  std::cout << "\n=== Final parameters (after command-line overrides) ===" << std::endl;
  std::cout << "  numberOfEvents : " << getParamS("numberOfEvents") << std::endl;
  std::cout << "  eventsInNtuple : " << getParamS("eventsInNtuple") << std::endl;
  std::cout << "  jobName        : " << getParamS("jobName")        << std::endl;
  std::cout << "  runNumber      : " << getParamS("runNumber")      << std::endl;
  std::cout << "  gun_particle   : " << getParamS("gun_particle",  false, "(default)") << std::endl;
  std::cout << "  gun_energy_min : " << getParamS("gun_energy_min",false, "(default)") << std::endl;
  std::cout << "  gun_energy_max : " << getParamS("gun_energy_max",false, "(default)") << std::endl;
  std::cout << "  sipmType       : " << getParamS("sipmType",      false, "(default)") << std::endl;
  std::cout << "======================================================\n" << std::endl;

  runConfig = getParamS("runConfig");
  runNumber = getParamI("runNumber");

  //  csv file defeinition.  (no CSV file in this program)
  // defineCSV("2dSC");
  // defineCSV("2dCH");
  // defineCSV("3dCH");

  //   root histogram/ntuple file definition
  string outname =
      getParamS("jobName") + "_run" + getParamS("runNumber") + "_" +
      getParamS("runSeq") + "_" + getParamS("runConfig") + "_" +
      getParamS("numberOfEvents") + "evt_";
  if (!getParamB("CaloXPythiaON", false, false))
  {
    outname += getParamS("gun_particle") + "_" + getParamS("gun_energy_min") + "_" + getParamS("gun_energy_max");
  }
  else
  {
    outname += "py8";
  }
  string outRootName = getParamS("rootPre") + "_" + outname + ".root";

  eventCounts = 0;
  eventCountsALL = 0;

  saveTruthHits = false;
  isMuon = false; // default is not muon
  if (getParamS("saveTruthHits").compare(0, 4, "true") == 0)
    saveTruthHits = true;

  //  ========  root histogram, ntuple file ===========
  fout = new TFile(outRootName.c_str(), "recreate");

  createNtuple = false;
  if (getParamS("createNtuple").compare(0, 4, "true") == 0)
    createNtuple = true;

  tree = new TTree("tree", "CaloX Tree");

  // set event counter.

  tree->Branch("run", &m_run);
  tree->Branch("event", &m_event);

  tree->Branch("beamMinE", &m_beamMinE);
  tree->Branch("beamMaxE", &m_beamMaxE);
  tree->Branch("gridSizeX", &m_gridSizeX);
  tree->Branch("gridSizeY", &m_gridSizeY);
  tree->Branch("gridSizeT", &m_gridSizeT);

  tree->Branch("calibSen", &m_calibSen);
  tree->Branch("calibSph", &m_calibSph);
  tree->Branch("calibCen", &m_calibCen);
  tree->Branch("calibCph", &m_calibCph);

  tree->Branch("beamX", &m_beamX);
  tree->Branch("beamY", &m_beamY);
  tree->Branch("beamZ", &m_beamZ);
  tree->Branch("beamE", &m_beamE);
  tree->Branch("beamID", &m_beamID);
  tree->Branch("beamType", &m_beamType);
  tree->Branch("beamPX", &m_beamPX);
  tree->Branch("beamPY", &m_beamPY);
  tree->Branch("beamPZ", &m_beamPZ);

  tree->Branch("ntruthhits", &m_nhitstruth);
  tree->Branch("truthhit_pid", &m_pidtruth);
  tree->Branch("truthhit_trackid", &m_trackidtruth);
  tree->Branch("truthhit_calotype", &m_calotypetruth);
  tree->Branch("truthhit_x", &m_xtruth);
  tree->Branch("truthhit_y", &m_ytruth);
  tree->Branch("truthhit_z", &m_ztruth);
  tree->Branch("truthhit_steplength", &m_steplengthtruth);
  tree->Branch("truthhit_globaltime", &m_globaltimetruth);
  tree->Branch("truthhit_localtime", &m_localtimetruth);
  tree->Branch("truthhit_edep", &m_edeptruth);
  tree->Branch("truthhit_edepNonIon", &m_edepNonIontruth);
  tree->Branch("truthhit_edepInv", &m_edepInvtruth);
  tree->Branch("truthhit_edepbirk", &m_edepbirktruth);
  tree->Branch("truthhit_process", &m_processtruth);
  tree->Branch("truthhit_ncer", &m_ncertruth);
  tree->Branch("truthhit_ncercap", &m_ncercaptruth);
  tree->Branch("truthhit_layerNumber", &m_layerNumber);
  tree->Branch("truthhit_rodNumber", &m_rodNumber);
  tree->Branch("truthhit_fiberNumber", &m_fiberNumber);

  tree->Branch("eCalotruth", &m_eCalotruth);
  tree->Branch("eWorldtruth", &m_eWorldtruth);
  tree->Branch("eLeaktruth", &m_eLeaktruth);
  tree->Branch("eLeaktruth_lt_z", &m_eLeaktruth_lt_z);
  tree->Branch("eLeaktruth_gt_z", &m_eLeaktruth_gt_z);
  tree->Branch("eInvisible", &m_eInvisible);
  tree->Branch("eRodtruth", &m_eRodtruth);
  tree->Branch("ePlatruth", &m_ePlatruth);
  tree->Branch("eQuatruth", &m_eQuatruth);
  tree->Branch("eScintruth", &m_eScintruth);

  tree->Branch("nhits3dSS", &m_nhits3dSS);
  tree->Branch("id3dSS", &m_id3dSS);
  tree->Branch("type3dSS", &m_type3dSS);
  tree->Branch("area3dSS", &m_area3dSS);
  tree->Branch("ix3dSS", &m_ix3dSS);
  tree->Branch("iy3dSS", &m_iy3dSS);
  tree->Branch("ixx3dSS", &m_ixx3dSS);
  tree->Branch("iyy3dSS", &m_iyy3dSS);
  tree->Branch("zslice3dSS", &m_zslice3dSS);
  tree->Branch("tslice3dSS", &m_tslice3dSS);
  tree->Branch("ph3dSS", &m_ph3dSS);
  tree->Branch("sum3dSS", &m_sum3dSS);

  tree->Branch("nhits3dCC", &m_nhits3dCC);
  tree->Branch("id3dCC", &m_id3dCC);
  tree->Branch("type3dCC", &m_type3dCC);
  tree->Branch("area3dCC", &m_area3dCC);
  tree->Branch("ix3dCC", &m_ix3dCC);
  tree->Branch("iy3dCC", &m_iy3dCC);
  tree->Branch("ixx3dCC", &m_ixx3dCC);
  tree->Branch("iyy3dCC", &m_iyy3dCC);
  tree->Branch("zslice3dCC", &m_zslice3dCC);
  tree->Branch("tslice3dCC", &m_tslice3dCC);
  tree->Branch("ph3dCC", &m_ph3dCC);
  tree->Branch("sum3dCC", &m_sum3dCC);

  tree->Branch("nhits3dQQ", &m_nhits3dQQ);
  tree->Branch("id3dQQ", &m_id3dQQ);
  tree->Branch("type3dQQ", &m_type3dQQ);
  tree->Branch("area3dQQ", &m_area3dQQ);
  tree->Branch("ix3dQQ", &m_ix3dQQ);
  tree->Branch("iy3dQQ", &m_iy3dQQ);
  tree->Branch("ixx3dQQ", &m_ixx3dQQ);
  tree->Branch("iyy3dQQ", &m_iyy3dQQ);
  tree->Branch("zslice3dQQ", &m_zslice3dQQ);
  tree->Branch("tslice3dQQ", &m_tslice3dQQ);
  tree->Branch("ph3dQQ", &m_ph3dQQ);
  tree->Branch("sum3dQQ", &m_sum3dQQ);

  tree->Branch("nOPs", &mP_nOPs);
  tree->Branch("OP_trackid", &mP_trackid);
  tree->Branch("OP_pos_produced_x", &mP_pos_produced_x);
  tree->Branch("OP_pos_produced_y", &mP_pos_produced_y);
  tree->Branch("OP_pos_produced_z", &mP_pos_produced_z);
  tree->Branch("OP_mom_produced_x", &mP_mom_produced_x);
  tree->Branch("OP_mom_produced_y", &mP_mom_produced_y);
  tree->Branch("OP_mom_produced_z", &mP_mom_produced_z);
  tree->Branch("OP_pos_final_x", &mP_pos_final_x);
  tree->Branch("OP_pos_final_y", &mP_pos_final_y);
  tree->Branch("OP_pos_final_z", &mP_pos_final_z);
  tree->Branch("OP_mom_final_x", &mP_mom_final_x);
  tree->Branch("OP_mom_final_y", &mP_mom_final_y);
  tree->Branch("OP_mom_final_z", &mP_mom_final_z);
  tree->Branch("OP_time_produced", &mP_time_produced);
  tree->Branch("OP_time_final", &mP_time_final);
  tree->Branch("OP_isCerenkov", &mP_isCerenkov);
  tree->Branch("OP_isScintillation", &mP_isScintillation);
  tree->Branch("OP_productionFiber", &mP_productionFiber);
  tree->Branch("OP_productionRod",   &mP_productionRod);
  tree->Branch("OP_productionLayer", &mP_productionLayer);
  tree->Branch("OP_finalFiber", &mP_finalFiber);
  tree->Branch("OP_isCoreC", &mP_isCoreC);
  tree->Branch("OP_isCoreS", &mP_isCoreS);
  tree->Branch("OP_isCoreQ", &mP_isCoreQ);
  tree->Branch("OP_isCladC", &mP_isCladC);
  tree->Branch("OP_isCladS", &mP_isCladS);
  tree->Branch("OP_isCladQ", &mP_isCladQ);
  tree->Branch("OP_pol_x", &mP_pol_x);
  tree->Branch("OP_pol_y", &mP_pol_y);
  tree->Branch("OP_pol_z", &mP_pol_z);

  tree->Branch("nOPsCer",     &mP_nOPsCer);
  tree->Branch("nOPsCer_Pla", &mP_nOPsCer_Pla);
  tree->Branch("nOPsCer_Qua", &mP_nOPsCer_Qua);
  tree->Branch("nOPsCer_Sci", &mP_nOPsCer_Sci);

  // Meridional (pz-only) analytical result
  tree->Branch("OP_isCaptured_m",            &mP_isCaptured_m);
  tree->Branch("OP_isAttenuated_m",          &mP_isAttenuated_m);
  tree->Branch("OP_captureAngle_m",          &mP_captureAngle_m);
  tree->Branch("OP_analyticalArrivalTime_m", &mP_analyticalArrivalTime_m);
  // Skew-corrected (L_z) analytical result
  tree->Branch("OP_isCaptured_s",            &mP_isCaptured_s);
  tree->Branch("OP_isAttenuated_s",          &mP_isAttenuated_s);
  tree->Branch("OP_captureAngle_s",          &mP_captureAngle_s);
  tree->Branch("OP_analyticalArrivalTime_s", &mP_analyticalArrivalTime_s);
}

// ########################################################################
CaloXTree::~CaloXTree() { std::cout << "deleting CaloXTree..." << std::endl; }

// ########################################################################
void CaloXTree::BeginEvent()
{
  // This is clled from PrimaryGeneratorAction::GeneratePrimaries,
  // not from CaloXEventAction..
  // clearCaloXHits();
  clearCaloXTree();
}

// ########################################################################
void CaloXTree::EndEvent()
{

  // std::cout<<"CaloXTree::EndEvent()  starting..."<<std::endl;
  eventCountsALL = eventCountsALL + 1;
  eventCounts = eventCounts + 1;
  mEvent = eventCounts;

  m_run = 1;
  m_event = eventCounts;

  if ((eventCounts - 1) < getParamI("eventsInNtuple"))
  {
    m_beamMinE = getParamF("gun_energy_min", false, -1.0);
    m_beamMaxE = getParamF("gun_energy_max", false, -1.0);
    m_gridSizeX = getParamF("gridSizeX"); // rod count
    m_gridSizeY = getParamF("gridSizeY"); // rod count
    m_gridSizeT = getParamF("gridSizeT"); // ps unit

    m_caloRotationX = getParamF("caloRotationX");
    m_caloRotationY = getParamF("caloRotationY");

    m_calibSen = getParamF("calibSen");
    m_calibSph = getParamF("calibSph");
    m_calibCen = getParamF("calibCen");
    m_calibCph = getParamF("calibCph");

    m_beamX = beamX;
    m_beamY = beamY;
    m_beamZ = beamZ;
    m_beamPX = beamPX;
    m_beamPY = beamPY;
    m_beamPZ = beamPZ;
    m_beamE = beamE;
    m_beamID = beamID;
    m_beamType = beamType;

    //  CC:  Cherenkov hits (ncer)
    m_sum3dCC = 0.0;
    for (auto itr = ctHits.begin(); itr != ctHits.end(); itr++)
    {
      CaloXID id(itr->first);
      int area = id.area(); // 0=Al-block, 1=no-SiPM, 2=6mm, 3=3mm
      if (area < 2)
        continue;
      double ncer = itr->second;
      if (round(ncer) < 1.0)
        continue; // 1.0 cherenkov photon cut
      m_sum3dCC = m_sum3dCC + ncer;
      // m_ky3dCC.push_back(itr->first);  // this used for debugging.
      int ky = id.iy() * 10; // 6mm SiPM
      if (area == 3)
      {
        ky = ky + id.iyy() + 1;
      } // 3mm SiPM
      m_id3dCC.push_back(id.ix() * 10000000 + ky * 1000 + id.tslice());
      m_type3dCC.push_back(id.type());
      m_area3dCC.push_back(id.area());
      m_ix3dCC.push_back(id.ix());
      m_iy3dCC.push_back(id.iy());
      m_ixx3dCC.push_back(id.ixx());
      m_iyy3dCC.push_back(id.iyy());
      m_zslice3dCC.push_back(id.zslice());
      m_tslice3dCC.push_back(id.tslice());
      m_ph3dCC.push_back(round(ncer));
    }
    m_nhits3dCC = m_ph3dCC.size();

    //  QQ: Quartz Cherenkov hits (ncer)
    m_sum3dQQ = 0.0;
    for (auto itr = qtHits.begin(); itr != qtHits.end(); itr++)
    {
      CaloXID id(itr->first);
      int area = id.area();
      if (area < 2)
        continue;
      double ncer = itr->second;
      if (round(ncer) < 1.0)
        continue;
      m_sum3dQQ = m_sum3dQQ + ncer;
      int ky = id.iy() * 10;
      if (area == 3)
      {
        ky = ky + id.iyy() + 1;
      }
      m_id3dQQ.push_back(id.ix() * 10000000 + ky * 1000 + id.tslice());
      m_type3dQQ.push_back(id.type());
      m_area3dQQ.push_back(id.area());
      m_ix3dQQ.push_back(id.ix());
      m_iy3dQQ.push_back(id.iy());
      m_ixx3dQQ.push_back(id.ixx());
      m_iyy3dQQ.push_back(id.iyy());
      m_zslice3dQQ.push_back(id.zslice());
      m_tslice3dQQ.push_back(id.tslice());
      m_ph3dQQ.push_back(round(ncer));
    }
    m_nhits3dQQ = m_ph3dQQ.size();

    //  SS: Scintillation hits (edepbirk)...
    m_sum3dSS = 0.0;
    for (auto itr = stHits.begin(); itr != stHits.end(); itr++)
    {
      CaloXID id(itr->first);
      int area = id.area(); // 0=Al-block, 1=no-SiPM, 2=6mm, 3=3mm
      if (area < 2)
        continue;
      double edepbirk = itr->second;
      if (edepbirk < 0.0001)
        continue; // 0.1 kev cut
      m_sum3dSS = m_sum3dSS + edepbirk;
      int ky = id.iy() * 10; // 6mm SiPM
      if (area == 3)
      {
        ky = ky + id.iyy() + 1;
      } // 3mm SiPM
      m_id3dSS.push_back(id.ix() * 10000000 + ky * 1000 + id.tslice());
      m_type3dSS.push_back(id.type());
      m_area3dSS.push_back(id.area());
      m_ix3dSS.push_back(id.ix());
      m_iy3dSS.push_back(id.iy());
      m_ixx3dSS.push_back(id.ixx());
      m_iyy3dSS.push_back(id.iyy());
      m_zslice3dSS.push_back(id.zslice());
      m_tslice3dSS.push_back(id.tslice());
      m_ph3dSS.push_back(edepbirk);
    }
    m_nhits3dSS = m_ph3dSS.size();

    m_nhitstruth = m_pidtruth.size();

    // optical photon hits
    for (auto const photon : photonData)
    {
      // if (photon.exitTime == 0.0)
      //   continue;
      mP_trackid.push_back(photon.trackID);
      mP_pos_produced_x.push_back(photon.productionPosition.x());
      mP_pos_produced_y.push_back(photon.productionPosition.y());
      mP_pos_produced_z.push_back(photon.productionPosition.z());
      mP_mom_produced_x.push_back(photon.productionMomentum.x());
      mP_mom_produced_y.push_back(photon.productionMomentum.y());
      mP_mom_produced_z.push_back(photon.productionMomentum.z());
      mP_pos_final_x.push_back(photon.exitPosition.x());
      mP_pos_final_y.push_back(photon.exitPosition.y());
      mP_pos_final_z.push_back(photon.exitPosition.z());
      mP_mom_final_x.push_back(photon.exitMomentum.x());
      mP_mom_final_y.push_back(photon.exitMomentum.y());
      mP_mom_final_z.push_back(photon.exitMomentum.z());
      mP_time_produced.push_back(photon.productionTime);
      mP_time_final.push_back(photon.exitTime);
      mP_isCerenkov.push_back(photon.isCerenkov);
      mP_isScintillation.push_back(photon.isScintillation);
      mP_productionFiber.push_back(photon.productionFiber);
      mP_productionRod.push_back(photon.productionRod);
      mP_productionLayer.push_back(photon.productionLayer);
      mP_finalFiber.push_back(photon.exitFiber);
      mP_isCoreC.push_back(photon.isCoreC);
      mP_isCoreS.push_back(photon.isCoreS);
      mP_isCoreQ.push_back(photon.isCoreQ);
      mP_isCladC.push_back(photon.isCladC);
      mP_isCladS.push_back(photon.isCladS);
      mP_isCladQ.push_back(photon.isCladQ);
      mP_pol_x.push_back(photon.polarization.x());
      mP_pol_y.push_back(photon.polarization.y());
      mP_pol_z.push_back(photon.polarization.z());
      mP_isCaptured_m.push_back(photon.isCaptured_m ? 1 : 0);
      mP_isAttenuated_m.push_back(photon.isAttenuated_m ? 1 : 0);
      mP_captureAngle_m.push_back(photon.captureAngle_m);
      mP_analyticalArrivalTime_m.push_back(photon.analyticalArrivalTime_m);
      mP_isCaptured_s.push_back(photon.isCaptured_s ? 1 : 0);
      mP_isAttenuated_s.push_back(photon.isAttenuated_s ? 1 : 0);
      mP_captureAngle_s.push_back(photon.captureAngle_s);
      mP_analyticalArrivalTime_s.push_back(photon.analyticalArrivalTime_s);

      // std::cout << "Propagation length in z " << photon.exitPosition.z() - photon.productionPosition.z() << " speed " << (photon.exitTime - photon.productionTime) / (photon.exitPosition.z() - photon.productionPosition.z()) << " costheta " << photon.productionMomentum.z() / photon.productionMomentum.mag() << std::endl;
    }
    mP_nOPs = photonData.size();

    //
    tree->Fill();
    std::cout << "Look into energy deposition in the calorimeter..." << std::endl;
    std::cout << "  eCalo=" << m_eCalotruth << "  eWorld=" << m_eWorldtruth << "  eLeak=" << m_eLeaktruth << "  eInvisible=" << m_eInvisible << "  eRod=" << m_eRodtruth << "  ePla(Plastic)=" << m_ePlatruth << "  eQua(Quartz)=" << m_eQuatruth << "  eScin=" << m_eScintruth << " eCalo+eWorld+eLeak+eInvisible=" << (m_eCalotruth + m_eWorldtruth + m_eLeaktruth + m_eInvisible) << std::endl;
  } //  end of if((eventCounts-1)<getParamI("eventsInNtuple"))
}

// ########################################################################
void CaloXTree::EndJob()
{
  fout->Write();
  fout->Close();
}
// ########################################################################
void CaloXTree::saveBeamXYZEPxPyPz(string ptype, int pdgid, float x, float y, float z,
                            float en, float px, float py, float pz)
{
  beamType = ptype; // sting pi+. e+ mu+ etc.
  beamID = pdgid;
  beamX = x; // in mm
  beamY = y;
  beamZ = z;
  beamE = en;
  beamPX = px;
  beamPY = py;
  beamPZ = pz;

  isMuon = (abs(pdgid) == 13);
}

// ########################################################################
void CaloXTree::clearCaloXTree()
{
  rtHits.clear();
  stHits.clear();
  ctHits.clear();
  qtHits.clear();
  rzHits.clear();
  szHits.clear();
  czHits.clear();
  qzHits.clear();
  rzEdep.clear();
  szEdep.clear();
  czEdep.clear();
  qzEdep.clear();
  // mEventNumber.clear();
  // mNxCell=0;
  // mNyCell=0;
  // mHepPy.clear();     // GeV

  m_beamX = 0.0;
  m_beamY = 0.0;
  m_beamZ = 0.0;
  m_beamPX = 0.0;
  m_beamPY = 0.0;
  m_beamPZ = 0.0;
  m_beamE = 0.0;
  m_beamID = 0;
  m_beamType = " ";
  isMuon = false;

  m_nhitstruth = 0;
  m_pidtruth.clear();
  m_trackidtruth.clear();
  m_calotypetruth.clear();
  m_xtruth.clear();
  m_ytruth.clear();
  m_ztruth.clear();
  m_steplengthtruth.clear();
  m_globaltimetruth.clear();
  m_localtimetruth.clear();
  m_edeptruth.clear();
  m_edepNonIontruth.clear();
  m_edepInvtruth.clear();
  m_edepbirktruth.clear();
  m_processtruth.clear();
  m_ncertruth.clear();
  m_ncercaptruth.clear();
  m_layerNumber.clear();
  m_rodNumber.clear();
  m_fiberNumber.clear();

  m_eCalotruth = 0.0;
  m_eWorldtruth = 0.0;
  m_eLeaktruth = 0.0;
  m_eLeaktruth_lt_z = 0.0;
  m_eLeaktruth_gt_z = 0.0;
  m_eInvisible = 0.0;
  m_eRodtruth = 0.0;
  m_ePlatruth = 0.0;
  m_eQuatruth = 0.0;
  m_eScintruth = 0.0;

  m_nhits3dSS = 0;
  m_id3dSS.clear();
  m_type3dSS.clear();
  m_area3dSS.clear();
  m_ix3dSS.clear();
  m_iy3dSS.clear();
  m_ixx3dSS.clear();
  m_iyy3dSS.clear();
  m_zslice3dSS.clear();
  m_tslice3dSS.clear();
  m_ph3dSS.clear();
  m_sum3dSS = 0.0;

  m_nhits3dCC = 0;
  //  m_ky3dCC.clear();   // this used for debugging.
  m_id3dCC.clear();
  m_type3dCC.clear();
  m_area3dCC.clear();
  m_ix3dCC.clear();
  m_iy3dCC.clear();
  m_ixx3dCC.clear();
  m_iyy3dCC.clear();
  m_zslice3dCC.clear();
  m_tslice3dCC.clear();
  m_ph3dCC.clear();
  m_sum3dCC = 0.0;

  m_nhits3dQQ = 0;
  m_id3dQQ.clear();
  m_type3dQQ.clear();
  m_area3dQQ.clear();
  m_ix3dQQ.clear();
  m_iy3dQQ.clear();
  m_ixx3dQQ.clear();
  m_iyy3dQQ.clear();
  m_zslice3dQQ.clear();
  m_tslice3dQQ.clear();
  m_ph3dQQ.clear();
  m_sum3dQQ = 0.0;

  // clean photons
  photonData.clear();
  mP_nOPs = 0;
  mP_trackid.clear();
  mP_pos_produced_x.clear();
  mP_pos_produced_y.clear();
  mP_pos_produced_z.clear();
  mP_mom_produced_x.clear();
  mP_mom_produced_y.clear();
  mP_mom_produced_z.clear();
  mP_pos_final_x.clear();
  mP_pos_final_y.clear();
  mP_pos_final_z.clear();
  mP_mom_final_x.clear();
  mP_mom_final_y.clear();
  mP_mom_final_z.clear();
  mP_time_produced.clear();
  mP_time_final.clear();
  mP_isCerenkov.clear();
  mP_isScintillation.clear();
  mP_productionFiber.clear();
  mP_productionRod.clear();
  mP_productionLayer.clear();
  mP_finalFiber.clear();
  mP_isCoreC.clear();
  mP_isCoreS.clear();
  mP_isCoreQ.clear();
  mP_isCladC.clear();
  mP_isCladS.clear();
  mP_isCladQ.clear();
  mP_pol_x.clear();
  mP_pol_y.clear();
  mP_pol_z.clear();
  mP_isCaptured_m.clear();
  mP_isAttenuated_m.clear();
  mP_captureAngle_m.clear();
  mP_analyticalArrivalTime_m.clear();
  mP_isCaptured_s.clear();
  mP_isAttenuated_s.clear();
  mP_captureAngle_s.clear();
  mP_analyticalArrivalTime_s.clear();

  mP_nOPsCer = 0;
  mP_nOPsCer_Pla = 0;
  mP_nOPsCer_Qua = 0;
  mP_nOPsCer_Sci = 0;
}


// ########################################################################
void CaloXTree::accumulateHits(CaloXHit ah)
{
  if (saveTruthHits && ah.edep >= 1.0e-6 && (ah.calotype > 1 || (isMuon && ah.calotype == 1)))
  {
    // save the truth hit in the scintillating and cherenkov fibers.
    // larger than 1 eV
    m_pidtruth.push_back(ah.pid);
    m_trackidtruth.push_back(ah.trackid);
    m_calotypetruth.push_back(ah.calotype);
    m_xtruth.push_back(ah.x);
    m_ytruth.push_back(ah.y);
    m_ztruth.push_back(ah.z);
    m_steplengthtruth.push_back(ah.steplength);
    m_globaltimetruth.push_back(ah.globaltime);
    m_localtimetruth.push_back(ah.localtime);
    m_edeptruth.push_back(ah.edep);
    m_edepNonIontruth.push_back(ah.edepNonIon);
    m_edepInvtruth.push_back(ah.edepInv);
    m_edepbirktruth.push_back(ah.edepbirk);
    m_processtruth.push_back(ah.process);
    m_ncertruth.push_back(ah.ncer);
    m_ncercaptruth.push_back(ah.ncercap);
    m_layerNumber.push_back(ah.layerNumber);
    m_rodNumber.push_back(ah.rodNumber);
    m_fiberNumber.push_back(ah.fiberNumber);
  }

  CaloXID id = ah.caloid;

  //   ROD;
  if (id.type() == 1)
  {
    int t = id.getTkey();
    rtHits[t] = rtHits[t] + ah.edep;
    int z = id.getZkey();
    rzHits[z] = rzHits[z] + ah.edep;
    rzEdep[z] = rzEdep[z] + ah.edep;
  }

  //   S-Fibers
  if (id.type() == 2)
  {
    int t = id.getTkey();
    // int t=2;
    stHits[t] = stHits[t] + ah.edepbirk;
    int z = id.getZkey();
    // int z=12;
    szHits[z] = szHits[z] + ah.edepbirk;
    szEdep[z] = szEdep[z] + ah.edep;
  }

  //   Plastic C-Fibers
  if (id.type() == 3)
  {
    int t = id.getTkey();
    ctHits[t] = ctHits[t] + ah.ncercap;
    int z = id.getZkey();
    czHits[z] = czHits[z] + ah.ncercap;
    czEdep[z] = czEdep[z] + ah.edep;
  }

  //   Quartz C-Fibers
  if (id.type() == 4)
  {
    int t = id.getTkey();
    qtHits[t] = qtHits[t] + ah.ncercap;
    int z = id.getZkey();
    qzHits[z] = qzHits[z] + ah.ncercap;
    qzEdep[z] = qzEdep[z] + ah.edep;
  }

  // mEventNumber.clear();
  // mNxCell=0;
  // mNyCell=0;
  // mHepPy.clear();     // GeV
}

void CaloXTree::accumulateEnergy(double edep, int type = 0)
{
  if (type == -99)
    m_eLeaktruth += edep;
  if (type == -91)
    m_eLeaktruth_gt_z += edep;
  if (type == -92)
    m_eLeaktruth_lt_z += edep;
  if (type == -90)
    m_eInvisible += edep;
  if (type == -1)
    m_eWorldtruth += edep;
  if (type >= 0)
    m_eCalotruth += edep;
  if (type == 1)
    m_eRodtruth += edep;
  if (type == 2)
    m_eScintruth += edep;
  if (type == 3)
    m_ePlatruth += edep;
  if (type == 4)
    m_eQuatruth += edep;
}

void CaloXTree::accumulateOPsCer(int fiberType, int nOPs)
{
  // fiberType: 0=scintillation, 1=plastic Cherenkov, 2=quartz Cherenkov
  mP_nOPsCer += nOPs;
  if (fiberType == 1)
    mP_nOPsCer_Pla += nOPs;
  else if (fiberType == 2)
    mP_nOPsCer_Qua += nOPs;
  else
    mP_nOPsCer_Sci += nOPs;
}

// ########################################################################
void CaloXTree::analyze()
{
}

//  =============================================================================
map<int, double> CaloXTree::make2Dhits(map<int, double> hits)
{
  // this produces 2D hits from 3d hits
  map<int, double> hits2d;

  for (auto itr = hits.begin(); itr != hits.end(); itr++)
  {
    // std::cout<<" rt-key "<<itr->first<<"  val "<<itr->second<<std::endl;
    int k = itr->first;
    int mask =
        (1 << 22) - 1; // assuming bit pattern {2,2,5,5,3,3,2,8} (see CaloXID)
    int newkey = k & mask;
    double val = itr->second;

    hits2d[newkey] = hits2d[newkey] + val;
  }
  return hits2d;
}

//  =============================================================================
void CaloXTree::defineCSV(string type)
{
  //  type: "2dSC", "2dCH", "3dCH"
  string csvname =
      "hit_" + type + "_" + getParamS("runNumber") + "_" + getParamS("runSeq") +
      "_" + getParamS("runConfig") + "_" + getParamS("numberOfEvents") + "ev_" +
      getParamS("gun_particle") + "_" + getParamS("gun_energy_min") + "_" +
      getParamS("gun_energy_max") + "_csv" + getParamS("csvHits" + type) +
      ".csv";

  //   open a csv file...
  fcsv[type] = std::make_unique<std::ofstream>(csvname, std::ios::out);

  if (fcsv[type]->is_open())
  {

    *(fcsv[type]) << "a,gun_particle,gun_energy_min,gun_energy_max"
                  << ",gridSizeX,gridSizeY,calibSen,calibCen,calibCph\n"
                  << "b"
                  << "," << getParamS("gun_particle") << ","
                  << getParamS("gun_energy_min") << ","
                  << getParamS("gun_energy_max") << ","
                  << getParamI("gridSizeX") << "," << getParamI("gridSizeY")
                  << "," << getParamF("calibSen") << ","
                  << getParamF("calibCen") << "," << getParamF("calibCph")
                  << "\n"
                  << "c"
                  << ",event"
                  << ",x"
                  << ",y"
                  << ",z"
                  << ",Ebeam"
                  << ",pid"
                  << ",totalphotons"
                  << ",npoints"
                  << ",ix"
                  << ",iy"
                  << ",photons"
                  << ",repeat(chID,photons)"
                  << "\n";
  }
  // }
}

//  =============================================================================
void CaloXTree::writeCSV(string type, map<int, double> &hits)
{
  //  type: "2dSC", "2dCH", "3dCH"

  // calculate sum and count hits...
  double sum = 0.0;
  int nhits = 0;
  for (auto itr = hits.begin(); itr != hits.end(); itr++)
  {
    double val = itr->second;
    sum = sum + val;
    nhits++;
  }

  // cout<<"CaloXTree::writeCSV  type "<<type<<endl;
  // for(auto itr=fcsv.begin(); itr !=fcsv.end(); itr++) {
  //    std::cout<<" irt-key "<<itr->first<<std::endl;
  // }

  if (!(fcsv[type]->is_open()))
  {
    cout << "csv file " << type << " is not open.  return without writing."
         << endl;
  }

  *(fcsv[type]) << eventCountsALL << "," << beamX << "," << beamY << ","
                << beamZ << "," << beamE << "," << beamType << "," << sum << ","
                << nhits;

  for (const auto &n : hits)
  { // using C__11 definition
    int key = n.first;
    CaloXID id = CaloXID(key);
    int ix = id.ix();
    int iy = id.iy();
    int it = id.tslice();
    int nph = int(n.second);
    *(fcsv[type]) << "," << ix << "," << iy << "," << it << "," << nph;
  }
  *(fcsv[type]) << endl;
}

//  =============================================================================
//  ============================================================================
//  === user run time parameter handling
//  =============================================================================
//  =============================================================================
void CaloXTree::readMacFile(string fileName)
{
  // read parameters from G4 mac file...
  ifstream macfile(fileName);

  if (macfile.is_open())
  {
    string line;
    while (getline(macfile, line))
    {
      // std::cout<<"line="<<line<<"=endline="<<std::endl;
      vector<string> tokens = parse_line(line);
      if (tokens.size() > 2)
      {
        for (unsigned i = 0; i < tokens.size(); i++)
        {
          if (tokens[0].compare("#$$$") == 0)
          {
            mcParams[tokens[1]] = tokens[2];
          }
        }
      }
    }
    macfile.close();
  }
  else
  {
    cout << "CaloXTree::readParamFile: error to open mac file, " << fileName
         << endl;
  }

  cout << " " << endl;
  cout << "=== Parameters from " << fileName
       << " (CaloXTree::readMacFile) ===" << endl;
  map<string, string>::iterator it;
  for (it = mcParams.begin(); it != mcParams.end(); ++it)
  {
    cout << it->first << " => " << it->second << endl;
  }
}

// =======================================================================
bool CaloXTree::setParam(string key, string val)
{
  bool keyfound = (mcParams.find(key) != mcParams.end());
  if (keyfound)
  {
    std::cout << "CaloXTree::setParam: (overwrite)  key=" << key
              << "  val=" << val << std::endl;
    mcParams[key] = val;
  }
  else
  {
    std::cout << "CaloXTree::setParam: key not found: " << key << " (inserting)" << std::endl;
    mcParams[key] = val;
  }
  return !keyfound; // return code: true means key was new (not an overwrite)
}

bool CaloXTree::getParamB(string key, bool is_required, bool default_value)
{
  // return true if the key is found in the map.
  // return false if the key is not found in the map.
  bool val = false;
  if (mcParams.find(key) != mcParams.end())
  {
    if (mcParams[key] == "true")
      val = true;
    else
      val = false;
  }
  else if (is_required)
  {
    std::cout << "  " << std::endl;
    std::cout << "CaloXTree::getParamB: Parameter key (" << key
              << ") does not exist in the mac file. Exit.." << std::endl;
    std::cout << "    note:  key word is case sensitive." << std::endl;
    std::cout << " bool " << std::endl;
    std::cout << "  " << std::endl;
    std::exit(0);
  }
  else
  {
    val = default_value; // return default value if not required.
  }
  return val;
}

// =======================================================================
float CaloXTree::getParamF(string key, bool is_required, float default_value)
{
  float val = 98765.0;
  if (mcParams.find(key) != mcParams.end())
  {
    val = std::stof(mcParams[key]);
  }
  else if (is_required)
  {
    std::cout << "  " << std::endl;
    std::cout << "CaloXTree::getParamF: Parameter key (" << key
              << ") does not exist in the mac file. Exit.." << std::endl;
    std::cout << "    note:  key word is case sensitive." << std::endl;
    std::cout << " float " << std::endl;
    std::cout << "  " << std::endl;
    std::exit(0);
  }
  else
  {
    val = default_value; // return default value if not required.
  }
  return val;
}

// =======================================================================
int CaloXTree::getParamI(string key, bool is_required, int default_value)
{
  int val = 98765;
  if (mcParams.find(key) != mcParams.end())
  {
    val = std::stoi(mcParams[key]);
  }
  else if (is_required)
  {
    std::cout << "  " << std::endl;
    std::cout << "CaloXTree::getParamI: Parameter key (" << key
              << ") does not exist in the mac file. Exit.." << std::endl;
    std::cout << "    note:  key word is case sensitive." << std::endl;
    std::cout << " int " << std::endl;
    std::cout << "  " << std::endl;
    std::exit(0);
  }
  else
  {
    val = default_value; // return default value if not required.
  }
  return val;
}

// =======================================================================
string CaloXTree::getParamS(string key, bool is_required, string default_value)
{
  string val = "aaa";
  if (mcParams.find(key) != mcParams.end())
  {
    val = mcParams[key];
  }
  else if (is_required)
  {
    std::cout << "  " << std::endl;
    std::cout << "CaloXTree::getParamS: Parameter key (" << key
              << ") does not exist in the mac file. Exit.." << std::endl;
    std::cout << "    note:  key word is case sensitive." << std::endl;
    std::cout << " string " << std::endl;
    std::cout << "  " << std::endl;
    std::exit(0);
  }
  else
  {
    val = default_value; // return default value if not required.
  }
  return val;
}

// =======================================================================
vector<string> CaloXTree::parse_line(string line)
{

  string buf;            // Have a buffer string
  stringstream ss(line); // Insert the string into a stream
  vector<string> tokens; // Create vector to hold our words

  if (line.size() > 0)
  {
    while (ss >> buf)
      tokens.push_back(buf);
  }
  return tokens;
  ;
};
