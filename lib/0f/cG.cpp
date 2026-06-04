#include "OperatorImport.h"
#include "complex_math.h"
#include "MSSM.h"

std::complex<double> MSSM::cG() {
    return (-0.16666666666666666*(pow(g3,3)*hbar)*EinsSum({LoopFunc({m3},8,mubarsq)},{{}},{}) + (2*pow(g3,3)*hbar)/5*EinsSum({LoopFunc({m3},12,mubarsq)},{{}},{}) - 0.05555555555555555*(pow(g3,3)*hbar)*EinsSum({LoopFunc({mdt},8,mubarsq)},{{1}},{}) + (pow(g3,3)*hbar)/8*EinsSum({LoopFunc({mdt},10,mubarsq)},{{1}},{}) - 0.06666666666666667*(pow(g3,3)*hbar)*EinsSum({LoopFunc({mdt},12,mubarsq)},{{1}},{}) - 0.1111111111111111*(pow(g3,3)*hbar)*EinsSum({LoopFunc({mqt},8,mubarsq)},{{1}},{}) + (pow(g3,3)*hbar)/4*EinsSum({LoopFunc({mqt},10,mubarsq)},{{1}},{}) + (-2*pow(g3,3)*hbar)/15*EinsSum({LoopFunc({mqt},12,mubarsq)},{{1}},{}) - 0.05555555555555555*(pow(g3,3)*hbar)*EinsSum({LoopFunc({mut},8,mubarsq)},{{1}},{}) + (pow(g3,3)*hbar)/8*EinsSum({LoopFunc({mut},10,mubarsq)},{{1}},{}) - 0.06666666666666667*(pow(g3,3)*hbar)*EinsSum({LoopFunc({mut},12,mubarsq)},{{1}},{}));
}
