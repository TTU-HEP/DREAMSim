//
// ---------------------------------------------------------------
//
// G4UserTrackingAction.cc
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "CaloXTrackingAction.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"

#include "G4Track.hh"

/////////////////////////////////////////////////////////
CaloXTrackingAction::CaloXTrackingAction()
    /////////////////////////////////////////////////////////
    : fpTrackingManager(0)
{
  if (!(G4ParticleTable::GetParticleTable()->GetReadiness()))
  {
    G4String msg;
    msg = " You are instantiating G4UserTrackingAction BEFORE your\n";
    msg += "G4VUserPhysicsList is instantiated and assigned to G4RunManager.\n";
    msg += " Such an instantiation is prohibited since Geant4 version 8.0. To fix this problem,\n";
    msg += "please make sure that your main() instantiates G4VUserPhysicsList AND\n";
    msg += "set it to G4RunManager before instantiating other user action classes\n";
    msg += "such as G4UserTrackingAction.";
    G4Exception("G4UserTrackingAction::G4UserTrackingAction()",
                "Tracking0001", FatalException, msg);
  }
}

/////////////////////////////////////////////////////////
CaloXTrackingAction::~CaloXTrackingAction()
/////////////////////////////////////////////////////////
{
  ;
}

//
void CaloXTrackingAction::PreUserTrackingAction(const G4Track *track)
{
  /*
   const G4DynamicParticle* dynamicParticle= track->GetDynamicParticle();
   G4int pdgcode=dynamicParticle->GetPDGcode();
   // G4int absPdgCode=abs(pdgcode);
   const G4ThreeVector& vtx=track->GetVertexPosition();

   if(track->GetTrackID() == 1) {
     std::cout<<"begin of track "<<track->GetTrackID()
            <<"   pID "<<pdgcode
            <<"   vtx x "<<vtx.x()<<" y "<<vtx.y()<<" z "<<vtx.z()
            <<std::endl;
   }
  */
}

//
void CaloXTrackingAction::PostUserTrackingAction(const G4Track *track)
{
}

/////////////////////////////////////////////////////////
void CaloXTrackingAction::
    SetTrackingManagerPointer(G4TrackingManager *pValue)
/////////////////////////////////////////////////////////
{
  fpTrackingManager = pValue;
}
