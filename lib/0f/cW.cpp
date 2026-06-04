#include "OperatorImport.h"
#include "complex_math.h"
#include "MSSM.h"

std::complex<double> MSSM::cW() {
    return (-0.1111111111111111*(pow(g2,3)*hbar)*EinsSum({LoopFunc({m2},8,mubarsq)},{{}},{}) + (4*pow(g2,3)*hbar)/15*EinsSum({LoopFunc({m2},12,mubarsq)},{{}},{}) - 0.05555555555555555*(pow(g2,3)*hbar)*EinsSum({LoopFunc({mPhi},8,mubarsq)},{{}},{}) + (pow(g2,3)*hbar)/8*EinsSum({LoopFunc({mPhi},10,mubarsq)},{{}},{}) - 0.06666666666666667*(pow(g2,3)*hbar)*EinsSum({LoopFunc({mPhi},12,mubarsq)},{{}},{}) - 0.05555555555555555*(pow(g2,3)*hbar)*EinsSum({LoopFunc({muTilde},8,mubarsq)},{{}},{}) + (2*pow(g2,3)*hbar)/15*EinsSum({LoopFunc({muTilde},12,mubarsq)},{{}},{}) - 0.05555555555555555*(pow(g2,3)*hbar)*EinsSum({LoopFunc({mlt},8,mubarsq)},{{1}},{}) + (pow(g2,3)*hbar)/8*EinsSum({LoopFunc({mlt},10,mubarsq)},{{1}},{}) - 0.06666666666666667*(pow(g2,3)*hbar)*EinsSum({LoopFunc({mlt},12,mubarsq)},{{1}},{}) - 0.16666666666666666*(pow(g2,3)*hbar)*EinsSum({LoopFunc({mqt},8,mubarsq)},{{1}},{}) + (3*pow(g2,3)*hbar)/8*EinsSum({LoopFunc({mqt},10,mubarsq)},{{1}},{}) - 0.2*(pow(g2,3)*hbar)*EinsSum({LoopFunc({mqt},12,mubarsq)},{{1}},{}));
}
