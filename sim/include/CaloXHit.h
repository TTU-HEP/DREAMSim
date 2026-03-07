#ifndef CaloXHit_h
#define CaloXHit_h 1

#include <vector>

#include "CaloXID.h"
#include "G4String.hh"

class CaloXHit
{
public:
   CaloXHit();
   ~CaloXHit();

   void print();

   CaloXID caloid;
   int pid;
   int trackid;
   int calotype;
   double x;
   double y;
   double z;
   double steplength;
   double globaltime;
   double localtime;
   double edep;
   double edepNonIon;
   double edepInv;
   double edepbirk;
   G4String process;
   double ncer;    // number of cerenkov photons
   double ncercap; // number of cerenkov photons
   int layerNumber;
   int rodNumber;
   int fiberNumber;
};

#endif
