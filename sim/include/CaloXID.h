#ifndef CaloXID_h
#define CaloXID_h 1

class CaloXID
{
public:
   CaloXID();
   ~CaloXID();
   CaloXID(int a_type, int a_fiber, int a_layer, int a_rod, double a_z, double a_t);
   CaloXID(int _key);

   int getTkey(); // time-slice (50ps/slice)  based key
   int getZkey(); // z-slice (cm/slice)) based key

   int type() { return _type; }
   int area() { return _area; }
   int ix() { return _ix; }
   int iy() { return _iy; }
   int ixx() { return _ixx; } //  sub address (0) for all.
   int iyy() { return _iyy; } //  sub-address (0,1,2,3) for 3mm SiPM in area 3.
   int zslice() { return _zslice; }
   int tslice() { return _tslice; }

   void print();

private:
   int findArea();

   int packKey(int type, int area, int ix, int iy, int ixx, int iyy, int ztype, int iz);
   int unpackKey(unsigned int k);

   float z0;
   float t0;
   float dz;
   float dt;
   int _layer; // [1,80] vertical
   int _rod;   // [1,90] horizontal axis
   int _fiber; // c[1,5], s[1,3]

   // key   (30 bits total)
   //   _type (2 bits)   1=rod, 2=sc,  3=ch
   //   _area (2 bits)   0=Al-block, 1=no-SiPM, 2=6mm, 3=3mm
   //   _ix   (5 bits)   [0,29]
   //   _iy   (5 bits)   [0,19]
   //   _ixx  (3 bits)   [0]
   //   _iyy  (3 bits)   [0,7]
   //   _ztype  (2 bit)  1=zslice, 2=tslice, 3=2D
   //   _iz   (8 bits)   [0,255]

   int _type; // [1,3]   1=rod, 2=sc, 3=cer
   int _area; //

   int _nx;     // =3: number of rods (in x) per channel
   int _ny;     // =4: number of layers (in y) per channel
   int _ix;     // rod/nx    [0,29]
   int _iy;     // layer/ny  [0,19]
   int _ztype;  // 1=zslice, 2=tslice, 3=2D
   int _zslice; //  [0,255]
   int _tslice; //  [0,255]

   int _nxx; //  =1: numbe of subchannels in x (3mm SiPM)
   int _nyy; //  =8: number of subchannel in y (3mm SiPM)
   int _ixx; //  sub channel id [0]  in area-3 (xx)
   int _iyy; //  sub channel id [0,7]  in area-3 (yy)
};

#endif
