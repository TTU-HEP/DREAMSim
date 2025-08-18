#ifndef CaloTree_h
#define CaloTree_h 1

#include <cstdlib> // for rand() on archer.
#include <fstream> // for input/output files
#include <iomanip> // for setw() in cout,
#include <iostream>
#include <map>
#include <math.h> // for sin(x) etc.
#include <memory>
#include <sstream> // for string stream
#include <string>
#include <vector>

class TFile;
class TTree;
class TH1D;
class TH2D;

class CaloHit;
struct PhotonInfo;

using namespace std;

class CaloTree
{
public:
  CaloTree(string, int argc, char **argv); // string outname
  ~CaloTree();                             // string outname
  void BeginEvent();
  void EndEvent();
  void EndJob();

  //  inout paramter handling ....
  bool setParam(string key, string val);
  bool getParamB(string key, bool is_required = true, bool default_value = false); // return true if the key is found in the map.
  float getParamF(string key, bool is_required = true, float default_value = 98765.0);
  double getParamD(string key, bool is_required = true, double default_value = 98765.0)
  {
    return static_cast<double>(getParamF(key, is_required, default_value));
  }
  int getParamI(string key, bool is_required = true, int default_value = 98765);
  string getParamS(string key, bool is_required = true, string default_value = "aaa");

  //  called fro SteppingAction...
  void accumulateHits(CaloHit aHit);
  void accumulateEnergy(double eleak, int type);
  void accumulateOPsCer(bool isCoreC = 0, int nOPs = 1);
  void saveBeamXYZEPxPyPz(string, int, float, float, float, float, float, float, float);

  // for histogrming...
  std::string title;
  std::map<std::string, TH1D *> histo1D;
  std::map<std::string, TH1D *>::iterator histo1Diter;
  std::map<std::string, TH2D *> histo2D;
  std::map<std::string, TH2D *>::iterator histo2Diter;

  vector<PhotonInfo> photonData;

private:
  // private functions.
  void readMacFile(string);
  vector<string> parse_line(string line);
  // std::map<std::string, std::string> mcParams;  //  MC run time parameters.
  map<string, string> mcParams; //  MC run time parameters.

  //
  void clearCaloTree();
  void analyze();

  map<int, double> make2Dhits(map<int, double> hits);
  void defineCSV(string type);
  void writeCSV(string type, map<int, double> &hits);

  string beamType;
  int beamID;
  float beamX, beamY, beamZ, beamE, beamPX, beamPY, beamPZ;
  //
  string runConfig;
  int runNumber;

  int eventCounts;
  int eventCountsALL;

  bool createNtuple;

  bool saveTruthHits;
  bool isMuon;

  // hit data in csv file
  map<string, int> csvEvents; // number of events to be written to csv file.
  map<std::string, std::unique_ptr<std::ofstream>> fcsv;
  // map<string,*ofstream> csvfiles;

  // ofstream hit2DFile;  //
  // ofstream hit3DFile;  //

  // ntuple file definition...
  TFile *fout;
  TTree *tree;

  //  accumulated energyr of photons
  //  in rods
  map<int, double> rtHits; // T-slice  (nominal 50 ps/slicen), edep
  map<int, double> rzHits; // Z-slice  (nominal 2 cm/slice)  , edep
  map<int, double> rzEdep; // Z-slice  (nominal 2 cm/slice)  , edep

  // in scit fibers
  map<int, double> stHits; // T-slice  (nominal 50 ps/slicen), edep-birk
  map<int, double> szHits; // Z-slice  (nominal 2 cm/slice)  , edep-birk
  map<int, double> szEdep; // Z-slice  (nominal 2 cm/slice)  , edep

  // in sherenkov fibers
  map<int, double> ctHits; // T-slice  (nominal 50 ps/slicen), n-photons
  map<int, double> czHits; // Z-slice  (nominal 2 cm/slice)  , n-photons
  map<int, double> czEdep; // Z-slice  (nominal 2 cm/slice)  , edep

  int mRun;
  int mEvent;

  int m_run;
  int m_event;

  float m_beamMinE;  // GeV
  float m_beamMaxE;  // GeV
  float m_gridSizeX; // mm
  float m_gridSizeY; // mm
  float m_gridSizeT; // ps
  float m_caloRotationX;
  float m_caloRotationY;

  float m_calibCore;  // ph to GeV
  float m_calibHelix; // ph to GeV

  float m_calibSen;
  float m_calibSph;
  float m_calibCen;
  float m_calibCph;

  int m_beamID; //  pdg ID
  float m_beamX;
  float m_beamY;
  float m_beamZ;
  float m_beamE;
  float m_beamPX; // GeV
  float m_beamPY; // GeV
  float m_beamPZ; // GeV
  string m_beamType;

  // truth hit variables (no sipm or time or position smearing applied)
  int m_nhitstruth;
  vector<int> m_pidtruth;
  vector<int> m_trackidtruth;
  vector<int> m_calotypetruth;
  vector<double> m_xtruth;
  vector<double> m_ytruth;
  vector<double> m_ztruth;
  vector<double> m_steplengthtruth;
  vector<double> m_globaltimetruth;
  vector<double> m_localtimetruth;
  vector<double> m_edeptruth;
  vector<double> m_edepNonIontruth;
  vector<double> m_edepInvtruth;
  vector<double> m_edepbirktruth;
  vector<std::string> m_processtruth;
  vector<double> m_ncertruth;
  vector<double> m_ncercaptruth;
  vector<int> m_layerNumber;
  vector<int> m_rodNumber;
  vector<int> m_fiberNumber;

  double m_eCalotruth;
  double m_eWorldtruth;
  double m_eLeaktruth;
  double m_eLeaktruth_gt_z;
  double m_eLeaktruth_lt_z;
  double m_eInvisible;
  double m_eRodtruth;
  double m_eCentruth;
  double m_eScintruth;

  // scintillation hit variables
  int m_nhits3dSS;
  vector<int> m_id3dSS; // channel ID   xxxyyyzzz
  vector<int> m_type3dSS;
  vector<int> m_area3dSS;
  vector<int> m_ix3dSS;
  vector<int> m_iy3dSS;
  vector<int> m_ixx3dSS;
  vector<int> m_iyy3dSS;
  vector<int> m_zslice3dSS;
  vector<int> m_tslice3dSS;
  vector<float> m_ph3dSS; // number of photons
  float m_sum3dSS;

  // cherenkov hit variables
  int m_nhits3dCC;
  vector<int> m_id3dCC; //  chnanel ID  xxxyyyttt
  vector<int> m_type3dCC;
  vector<int> m_area3dCC;
  vector<int> m_ix3dCC;
  vector<int> m_iy3dCC;
  vector<int> m_ixx3dCC;
  vector<int> m_iyy3dCC;
  vector<int> m_zslice3dCC;
  vector<int> m_tslice3dCC;
  vector<float> m_ph3dCC; //  number of photons
  float m_sum3dCC;

  // optical photon hit variables
  int mP_nOPs;
  vector<int> mP_trackid;
  vector<double> mP_pos_produced_x;
  vector<double> mP_pos_produced_y;
  vector<double> mP_pos_produced_z;
  vector<double> mP_mom_produced_x;
  vector<double> mP_mom_produced_y;
  vector<double> mP_mom_produced_z;
  vector<double> mP_time_produced;
  vector<double> mP_pos_final_x;
  vector<double> mP_pos_final_y;
  vector<double> mP_pos_final_z;
  vector<double> mP_mom_final_x;
  vector<double> mP_mom_final_y;
  vector<double> mP_mom_final_z;
  vector<double> mP_time_final;
  vector<int> mP_productionFiber;
  vector<int> mP_finalFiber;
  vector<int> mP_isCerenkov;
  vector<int> mP_isScintillation;
  vector<bool> mP_isCoreC;
  vector<bool> mP_isCoreS;
  vector<bool> mP_isCladC;
  vector<bool> mP_isCladS;
  vector<double> mP_pol_x;
  vector<double> mP_pol_y;
  vector<double> mP_pol_z;

  int mP_nOPsCer;     // number of Cerenkov photons
  int mP_nOPsCer_Cer; // number of Cerenkov photons in cherenkov fibers
  int mP_nOPsCer_Sci; // number of Cerenkov photons in scintillation fibers
};

#endif
