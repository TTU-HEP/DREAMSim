#include "CaloXFiberMap.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

namespace
{
   //  A deliberately small JSON reader.  The map file is written by
   //  tools/fibermap_from_image.py and only ever holds numbers, strings and
   //  arrays of strings, so pulling out the handful of keys we need by hand is
   //  simpler than carrying a JSON library around.  Keys are matched with their
   //  quotes so that "grid" does not also match "gridNote".

   size_t findKey(const std::string &s, const std::string &key, size_t from)
   {
      return s.find("\"" + key + "\"", from);
   }

   bool getInt(const std::string &s, const std::string &key, size_t from, int &out)
   {
      size_t k = findKey(s, key, from);
      if (k == std::string::npos)
         return false;
      size_t c = s.find(':', k);
      if (c == std::string::npos)
         return false;
      out = std::atoi(s.c_str() + c + 1);
      return true;
   }

   bool getStringArray(const std::string &s, const std::string &key, size_t from,
                       std::vector<std::string> &out)
   {
      size_t k = findKey(s, key, from);
      if (k == std::string::npos)
         return false;
      size_t open = s.find('[', k);
      size_t close = s.find(']', open);
      if (open == std::string::npos || close == std::string::npos)
         return false;
      out.clear();
      size_t p = open;
      while (true)
      {
         size_t q0 = s.find('"', p + 1);
         if (q0 == std::string::npos || q0 > close)
            break;
         size_t q1 = s.find('"', q0 + 1);
         if (q1 == std::string::npos || q1 > close)
            break;
         out.push_back(s.substr(q0 + 1, q1 - q0 - 1));
         p = q1;
      }
      return !out.empty();
   }
} // namespace

CaloXFiberMap::CaloXFiberMap()
    : fLoaded(false), fNRods(0), fNLayers(0), fRodsPerCellX(3), fLayersPerCellY(4),
      fIxMin(0), fIxMax(-1), fIyMin(0), fIyMax(-1)
{
}

CaloXFiberMap::Type CaloXFiberMap::fromChar(char c)
{
   if (c == 'Q' || c == 'q')
      return kQuartz;
   if (c == 'P' || c == 'p')
      return kPlastic;
   return kEmpty; //  '.', '_' and anything unexpected: no fibers
}

bool CaloXFiberMap::load(const std::string &fileName)
{
   fLoaded = false;

   std::ifstream in(fileName.c_str());
   if (!in)
   {
      std::cout << "CaloXFiberMap: cannot open fiber map file " << fileName << std::endl;
      return false;
   }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string s = ss.str();

   if (!getInt(s, "nRods", 0, fNRods) || !getInt(s, "nLayers", 0, fNLayers) ||
       !getInt(s, "rodsPerCellX", 0, fRodsPerCellX) ||
       !getInt(s, "layersPerCellY", 0, fLayersPerCellY))
   {
      std::cout << "CaloXFiberMap: " << fileName
                << " is missing one of nRods, nLayers, rodsPerCellX, layersPerCellY"
                << std::endl;
      return false;
   }
   if (fRodsPerCellX <= 0 || fLayersPerCellY <= 0)
   {
      std::cout << "CaloXFiberMap: rodsPerCellX and layersPerCellY must be positive"
                << std::endl;
      return false;
   }

   if (!getStringArray(s, "grid", 0, fGrid))
   {
      std::cout << "CaloXFiberMap: " << fileName << " has no grid" << std::endl;
      return false;
   }

   const int nCellsX = fNRods / fRodsPerCellX;
   const int nCellsY = fNLayers / fLayersPerCellY;
   if ((int)fGrid.size() != nCellsY)
   {
      std::cout << "CaloXFiberMap: grid has " << fGrid.size() << " rows but nLayers/"
                << fLayersPerCellY << " = " << nCellsY << std::endl;
      return false;
   }
   for (size_t r = 0; r < fGrid.size(); ++r)
   {
      if ((int)fGrid[r].size() != nCellsX)
      {
         std::cout << "CaloXFiberMap: grid row " << r << " has " << fGrid[r].size()
                   << " columns but nRods/" << fRodsPerCellX << " = " << nCellsX
                   << std::endl;
         return false;
      }
   }

   //  The central region is optional: without it the whole detector is read at
   //  the coarse 3 x 4 granularity.
   const size_t central = findKey(s, "central", 0);
   if (central != std::string::npos)
   {
      if (!getInt(s, "ixMin", central, fIxMin) || !getInt(s, "ixMax", central, fIxMax) ||
          !getInt(s, "iyMin", central, fIyMin) || !getInt(s, "iyMax", central, fIyMax) ||
          !getStringArray(s, "grid", central, fCentral))
      {
         std::cout << "CaloXFiberMap: the central block of " << fileName
                   << " is incomplete" << std::endl;
         return false;
      }
      const int wantRows = (fIyMax - fIyMin + 1) * fLayersPerCellY;
      const int wantCols = fIxMax - fIxMin + 1;
      if ((int)fCentral.size() != wantRows)
      {
         std::cout << "CaloXFiberMap: central grid has " << fCentral.size()
                   << " rows, expected " << wantRows << std::endl;
         return false;
      }
      for (size_t r = 0; r < fCentral.size(); ++r)
      {
         if ((int)fCentral[r].size() != wantCols)
         {
            std::cout << "CaloXFiberMap: central grid row " << r << " has "
                      << fCentral[r].size() << " columns, expected " << wantCols
                      << std::endl;
            return false;
         }
      }
   }

   fLoaded = true;
   return true;
}

CaloXFiberMap::Type CaloXFiberMap::type(int rod, int layer) const
{
   if (!fLoaded || rod < 0 || layer < 0 || rod >= fNRods || layer >= fNLayers)
      return kEmpty;

   const int ix = rod / fRodsPerCellX;
   const int iy = layer / fLayersPerCellY;
   const int nCellsY = fNLayers / fLayersPerCellY;

   //  Central cells are one layer tall, so they are indexed by the layer itself.
   if (!fCentral.empty() && ix >= fIxMin && ix <= fIxMax && iy >= fIyMin && iy <= fIyMax)
   {
      const int iyy = layer % fLayersPerCellY;
      const int row = (fIyMax - iy) * fLayersPerCellY + (fLayersPerCellY - 1 - iyy);
      const int col = ix - fIxMin;
      return fromChar(fCentral[row][col]);
   }

   return fromChar(fGrid[nCellsY - 1 - iy][ix]);
}

void CaloXFiberMap::print() const
{
   if (!fLoaded)
   {
      std::cout << "CaloXFiberMap: not loaded" << std::endl;
      return;
   }

   long nQ = 0, nP = 0, nE = 0;
   for (int layer = 0; layer < fNLayers; ++layer)
      for (int rod = 0; rod < fNRods; ++rod)
      {
         switch (type(rod, layer))
         {
         case kQuartz: ++nQ; break;
         case kPlastic: ++nP; break;
         default: ++nE; break;
         }
      }

   std::cout << "CaloXFiberMap: " << fNRods << " rods x " << fNLayers << " layers, cells of "
             << fRodsPerCellX << " x " << fLayersPerCellY << std::endl;
   if (!fCentral.empty())
      std::cout << "  central region (1 layer per cell): ix " << fIxMin << ".." << fIxMax
                << ", iy " << fIyMin << ".." << fIyMax << std::endl;
   std::cout << "  coppers: " << nQ << " quartz (4C+3S), " << nP << " plastic (4C+3S), "
             << nE << " empty" << std::endl;

   //  One character per cell, highest iy first, so the run log shows the map
   //  that was actually built.
   const int nCellsX = fNRods / fRodsPerCellX;
   const int nCellsY = fNLayers / fLayersPerCellY;
   for (int iy = nCellsY - 1; iy >= 0; --iy)
   {
      std::string line;
      for (int ix = 0; ix < nCellsX; ++ix)
      {
         const int rod = ix * fRodsPerCellX;
         bool mixed = false;
         Type t0 = type(rod, iy * fLayersPerCellY);
         for (int k = 1; k < fLayersPerCellY; ++k)
            if (type(rod, iy * fLayersPerCellY + k) != t0)
               mixed = true;
         if (mixed)
            line += '*'; //  central cells, which differ layer by layer
         else
            line += (t0 == kQuartz ? 'Q' : (t0 == kPlastic ? 'P' : '.'));
      }
      std::cout << "  " << line << "   iy=" << iy << std::endl;
   }
}
