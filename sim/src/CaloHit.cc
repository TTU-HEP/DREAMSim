#include "CaloHit.h"

#include <iostream> // for cout

CaloHit::CaloHit()
    : pid(0), trackid(0), calotype(0), x(0.0), y(0.0), z(0.0), steplength(0.0), globaltime(0.0), localtime(0.0), edep(0.0), edepNonIon(0.), edepInv(0.), edepbirk(0.0), process(""), ncer(0.0), ncercap(0.0), layerNumber(-1), rodNumber(-1), fiberNumber(-1)
{
}

CaloHit::~CaloHit()
{
}

// ------------------------------------------------------------------
void CaloHit::print()
{
    std::cout << "CaloHit pid " << pid << " xyz(" << x << "," << y << "," << z << ") L " << steplength;
    std::cout << " T " << globaltime << "  edep " << edep << "  ncer " << ncer;
    std::cout << " " << std::endl;
    caloid.print();
}
