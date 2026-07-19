/// \file CaloXRunAction.hh
/// \brief Definition of the CaloXRunAction class

#ifndef CaloXRunAction_h
#define CaloXRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Timer.hh"
#include "globals.hh"

class CaloXTree;
class CaloXRunAction : public G4UserRunAction
{
public:
  CaloXRunAction(CaloXTree *);
  virtual ~CaloXRunAction();

  virtual void BeginOfRunAction(const G4Run *);
  virtual void EndOfRunAction(const G4Run *);

private:
  CaloXTree *hh;
  G4Timer fTimer; // measures real (wall) and CPU time of the run
};


#endif
