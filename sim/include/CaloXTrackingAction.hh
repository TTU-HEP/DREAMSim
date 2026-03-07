//
//
//---------------------------------------------------------------
//
// G4UserTrackingAction.hh
//
// class description:
//   This class represents actions taken place by the user at
//   the start/end point of processing one track.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

class CaloXTrackingAction;

#ifndef CaloXTrackingAction_h
#define CaloXTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class G4TrackingManager; // Forward declaration
class G4Track;

///////////////////////////
class CaloXTrackingAction : public G4UserTrackingAction
///////////////////////////
{

   //--------
public: // with description
        //--------

   // Constructor & Destructor
   CaloXTrackingAction();
   virtual ~CaloXTrackingAction();

   // Member functions
   void SetTrackingManagerPointer(G4TrackingManager *pValue);
   void PreUserTrackingAction(const G4Track *);
   void PostUserTrackingAction(const G4Track *);

   //-----------
protected:
   //-----------

   // Member data
   G4TrackingManager *fpTrackingManager;
};

#endif
