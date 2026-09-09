#ifndef CaloXFiberMap_h
#define CaloXFiberMap_h 1

#include <string>
#include <vector>

//  Which fibers sit inside each copper rod.
//
//  The real detector is not uniform: a copper either carries 4 quartz Cherenkov
//  fibers plus 3 scintillating ones, or 4 plastic Cherenkov fibers plus 3
//  scintillating ones, or no fibers at all.  Which of the three it is, is given
//  per cell of 3 rods (x) by 4 layers (y) -- the same cells that CaloXID uses
//  for the readout, so cell (ix, iy) here is exactly CaloXID::ix(), iy().  In
//  the central region the cells are 3 rods by 1 layer instead, and are indexed
//  by (ix, iy, iyy) with iyy = layer % 4, again as in CaloXID.
//
//  The map is read from a JSON file written by tools/fibermap_from_image.py.
class CaloXFiberMap
{
public:
   enum Type
   {
      kEmpty = 0,  //  solid copper, no fibers
      kPlastic = 1, //  4 plastic Cherenkov + 3 scintillating
      kQuartz = 2  //  4 quartz Cherenkov  + 3 scintillating
   };

   CaloXFiberMap();

   //  Read the map.  Returns false (and leaves the map unloaded) on any problem,
   //  after printing what went wrong.
   bool load(const std::string &fileName);

   bool loaded() const { return fLoaded; }

   //  The fiber content of one copper.  Rods and layers outside the map, and
   //  any character the map does not recognise, come back as kEmpty.
   Type type(int rod, int layer) const;

   int nRods() const { return fNRods; }
   int nLayers() const { return fNLayers; }
   int rodsPerCellX() const { return fRodsPerCellX; }
   int layersPerCellY() const { return fLayersPerCellY; }

   //  Counts per type plus the map as it will be built, for the run log.
   void print() const;

private:
   bool fLoaded;
   int fNRods;
   int fNLayers;
   int fRodsPerCellX;
   int fLayersPerCellY;

   std::vector<std::string> fGrid;    //  one row per iy, highest iy first
   std::vector<std::string> fCentral; //  one row per layer, topmost layer first
   int fIxMin, fIxMax, fIyMin, fIyMax;

   static Type fromChar(char c);
};

#endif
