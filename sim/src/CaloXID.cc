#include "CaloXID.h"

#include <iostream> // for cout

CaloXID::CaloXID()
{
}

CaloXID::~CaloXID()
{
   // std::cout<<"sk deleting CaloXID..."<<std::endl;
}

CaloXID::CaloXID(int a_calotype, int a_fiber, int a_layer, int a_rod, double a_z, double a_t)
{
   _type = a_calotype;
   _layer = a_layer;
   _rod = a_rod;
   _fiber = a_fiber;

   _nx = 3; //   number of rods per channel
   _ny = 4; //   number of layers per channel
   _ix = _rod / _nx;
   _iy = _layer / _ny;

   _area = findArea();

   _ixx = 0;
   _iyy = 0;
   if (_area == 3)
   {
      _iyy = _layer % 4;
   }

   _zslice = 0;
   _tslice = 0;

   z0 = -1000.0; // mm
   dz = 20.0;    // mm
   _zslice = int((a_z - z0) / dz);
   if (_zslice < 0)
      _zslice = 0;
   if (_zslice > 500)
      _zslice = 500;

   t0 = 0.0;
   dt = 0.05;             // nsec
   double c = 300.0;      //  speed of light 300 mm/nsec
   double zback = 2000.0; //  length of the calorimeter
   // tslice=int((a_t-t0)/dt)+10;
   double zlocal = a_z - z0;
   double tlocal = a_t - t0;
   double tof = zlocal / c;
   double tback = (zback - zlocal) / (c * 19.0 / 30.0);
   // double t=(tlocal-tof)+tback;
   double t = tlocal + tback;
   _tslice = int(t / dt);
   if (_tslice < 0)
      _tslice = 0;
   if (_tslice > 511)
      _tslice = 511;
}

// ------------------------------------------------------------------------------------
int CaloXID::findArea()
{

   int a = 0;
   if (_ix > 4 && _ix < 25 && _iy > 1 && _iy < 18)
   {
      a = 2; // sipm 6 mm
   }
   if (_ix > 12 && _ix < 17)
   {
      if (_iy > 17 && _iy < 20)
      {
         a = 2; // sipm 6 mm
      }
      if (_iy > -1 && _iy < 3)
      {
         a = 2; // sipm 6 mm
      }
      if (_iy > 7 && _iy < 12)
      {
         a = 3; // sipm 6 mm
      }
   }

   return a;
}

// ------------------------------------------------------------------------------------
CaloXID::CaloXID(int _key)
{
   //  input _key is defined for only one case,  tslice or zslice...
   // dt=0.05 ;  // nsec
   // dz=20.0 ;  // mm
   // _type=(_key%10);
   // _layer=(_key%1000)/10;
   // _rod=(_key%100000)/1000;
   // int a=_key/100000;
   // if (a<500) {
   //    _iz=0;
   //   _it=a;
   // } else {
   //    _iz=a-500;
   //    _it=0;
   //}
   unpackKey(_key);
}

// ------------------------------------------------------------------------------------
int CaloXID::getTkey()
{
   // return _type + _iy*10 + _ix*1000 + _it*100000 ;
   return packKey(_type, _area, _ix, _iy, _ixx, _iyy, 2, _tslice);
}
// ------------------------------------------------------------------------------------
int CaloXID::getZkey()
{
   // return _type + _iy*10 + _ix*1000 + (_iz+500)*100000 ;
   return packKey(_type, _area, _ix, _iy, _ixx, _iyy, 1, _zslice);
}

// ------------------------------------------------------------------------------------
int CaloXID::packKey(int type, int area, int ix, int iy, int ixx, int iyy, int ztype, int iz)
{
   unsigned int sbits[8] = {2, 2, 5, 5, 3, 3, 2, 8}; // bits for packing (31 bits)

   unsigned int k = 0;
   unsigned int n = 0;
   k |= type << n;

   n = n + sbits[0];
   k |= area << n;

   n = n + sbits[1];
   k |= ix << n;

   n = n + sbits[2];
   k |= iy << n;

   n = n + sbits[3];
   k |= ixx << n;

   n = n + sbits[4];
   k |= iyy << n;

   n = n + sbits[5];
   k |= ztype << n;

   n = n + sbits[6];
   k |= iz << n;

   return k;
}
// ------------------------------------------------------------------------------------
int CaloXID::unpackKey(unsigned int k)
{

   unsigned int sbits[8] = {2, 2, 5, 5, 3, 3, 2, 8}; // bits for packing (31 bits)

   unsigned int n = 0;
   unsigned int mask = (1 << sbits[0]) - 1;
   _type = k & mask;

   n = n + sbits[0];
   mask = (1 << sbits[1]) - 1;
   _area = (k >> n) & mask;

   n = n + sbits[1];
   mask = (1 << sbits[2]) - 1;
   _ix = (k >> n) & mask;

   n = n + sbits[2];
   mask = (1 << sbits[3]) - 1;
   _iy = (k >> n) & mask;

   n = n + sbits[3];
   mask = (1 << sbits[4]) - 1;
   _ixx = (k >> n) & mask;

   n = n + sbits[4];
   mask = (1 << sbits[5]) - 1;
   _iyy = (k >> n) & mask;

   n = n + sbits[5];
   mask = (1 << sbits[6]) - 1;
   int ztype = (k >> n) & mask;

   n = n + sbits[6];
   mask = (1 << sbits[7]) - 1;
   int iz = (k >> n) & mask;

   if (ztype == 1)
   {
      _zslice = iz;
      _tslice = 0;
   }

   if (ztype == 2)
   {
      _zslice = 0;
      _tslice = iz;
   }

   if (ztype == 3)
   {
      _zslice = 0;
      _tslice = 0;
   }

   return 0;
}

// ------------------------------------------------------------------------------------
void CaloXID::print()
{
   std::cout << "CaloXID type " << _type << "  ix " << _ix << "  iy " << _iy;
   std::cout << "  sub(ixx " << _ixx << "  iyy " << _iyy << ")";
   std::cout << "  zslice " << _zslice << "  tslice " << _tslice << std::endl;
}
