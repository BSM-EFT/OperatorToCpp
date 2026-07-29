#include "LF.h"
#include "LF_helper.h"

// elaborate definitions of loop-functions in terms of masses
std::complex<double> LF(const std::vector<std::complex<double> >& m, int code, double mubarsq) {
    switch(code) {
        case 1: 
            return LFh1(m, mubarsq); 
            break;
        case 2: 
            return LFh2(m, mubarsq); 
            break;
        case 3: 
            return LFh3(m, mubarsq); 
            break;
        case 4: 
            return LFh4(m, mubarsq); 
            break;
        case 5: 
            return LFh5(m, mubarsq); 
            break;
        case 6: 
            return LFh6(m, mubarsq); 
            break;
        case 7: 
            return LFh7(m, mubarsq); 
            break;
        case 8: 
            return LFh8(m, mubarsq); 
            break;
        case 9: 
            return LFh9(m, mubarsq); 
            break;
        case 10: 
            return LFh10(m, mubarsq); 
            break;
        case 11: 
            return LFh11(m, mubarsq); 
            break;
        case 12: 
            return LFh12(m, mubarsq); 
            break;
        case 13: 
            return LFh13(m, mubarsq); 
            break;
        case 14: 
            return LFh14(m, mubarsq); 
            break;
        case 15: 
            return LFh15(m, mubarsq); 
            break;
        case 16: 
            return LFh16(m, mubarsq); 
            break;
        case 17: 
            return LFh17(m, mubarsq); 
            break;
        case 18: 
            return LFh18(m, mubarsq); 
            break;
        case 19: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},4,mubarsq); }
            else { return LFh19(m, mubarsq); }
            break;
        case 20: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},5,mubarsq); }
            else { return LFh20(m, mubarsq); }
            break;
        case 21: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},7,mubarsq); }
            else { return LFh21(m, mubarsq); }
            break;
        case 22: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},8,mubarsq); }
            else { return LFh22(m, mubarsq); }
            break;
        case 23: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},10,mubarsq); }
            else { return LFh23(m, mubarsq); }
            break;
        case 24: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},11,mubarsq); }
            else { return LFh24(m, mubarsq); }
            break;
        case 25: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},13,mubarsq); }
            else { return LFh25(m, mubarsq); }
            break;
        case 26: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},14,mubarsq); }
            else { return LFh26(m, mubarsq); }
            break;
        case 27: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},18,mubarsq); }
            else { return LFh27(m, mubarsq); }
            break;
        case 28: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},7,mubarsq); }
            else { return LFh28(m, mubarsq); }
            break;
        case 29: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},8,mubarsq); }
            else { return LFh29(m, mubarsq); }
            break;
        case 30: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},9,mubarsq); }
            else { return LFh30(m, mubarsq); }
            break;
        case 31: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},10,mubarsq); }
            else { return LFh31(m, mubarsq); }
            break;
        case 32: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},11,mubarsq); }
            else { return LFh32(m, mubarsq); }
            break;
        case 33: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},12,mubarsq); }
            else { return LFh33(m, mubarsq); }
            break;
        case 34: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},13,mubarsq); }
            else { return LFh34(m, mubarsq); }
            break;
        case 35: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},14,mubarsq); }
            else { return LFh35(m, mubarsq); }
            break;
        case 36: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},17,mubarsq); }
            else { return LFh36(m, mubarsq); }
            break;
        case 37: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},18,mubarsq); }
            else { return LFh37(m, mubarsq); }
            break;
        case 38: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},9,mubarsq); }
            else { return LFh38(m, mubarsq); }
            break;
        case 39: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},10,mubarsq); }
            else { return LFh39(m, mubarsq); }
            break;
        case 40: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},11,mubarsq); }
            else { return LFh40(m, mubarsq); }
            break;
        case 41: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},12,mubarsq); }
            else { return LFh41(m, mubarsq); }
            break;
        case 42: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},13,mubarsq); }
            else { return LFh42(m, mubarsq); }
            break;
        case 43: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},14,mubarsq); }
            else { return LFh43(m, mubarsq); }
            break;
        case 44: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},15,mubarsq); }
            else { return LFh44(m, mubarsq); }
            break;
        case 45: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},16,mubarsq); }
            else { return LFh45(m, mubarsq); }
            break;
        case 46: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},17,mubarsq); }
            else { return LFh46(m, mubarsq); }
            break;
        case 47: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},18,mubarsq); }
            else { return LFh47(m, mubarsq); }
            break;
        case 48: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},12,mubarsq); }
            else { return LFh48(m, mubarsq); }
            break;
        case 49: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},13,mubarsq); }
            else { return LFh49(m, mubarsq); }
            break;
        case 50: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},14,mubarsq); }
            else { return LFh50(m, mubarsq); }
            break;
        case 51: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},15,mubarsq); }
            else { return LFh51(m, mubarsq); }
            break;
        case 52: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},16,mubarsq); }
            else { return LFh52(m, mubarsq); }
            break;
        case 53: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},17,mubarsq); }
            else { return LFh53(m, mubarsq); }
            break;
        case 54: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},18,mubarsq); }
            else { return LFh54(m, mubarsq); }
            break;
        case 55: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},15,mubarsq); }
            else { return LFh55(m, mubarsq); }
            break;
        case 56: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},16,mubarsq); }
            else { return LFh56(m, mubarsq); }
            break;
        case 57: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},17,mubarsq); }
            else { return LFh57(m, mubarsq); }
            break;
        case 58: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0]},18,mubarsq); }
            else { return LFh58(m, mubarsq); }
            break;
        case 59: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},28,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},28,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},21,mubarsq); }
            else { return LFh59(m, mubarsq); }
            break;
        case 60: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},29,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},29,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},22,mubarsq); }
            else { return LFh60(m, mubarsq); }
            break;
        case 61: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},31,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},39,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},23,mubarsq); }
            else { return LFh61(m, mubarsq); }
            break;
        case 62: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},32,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},40,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},24,mubarsq); }
            else { return LFh62(m, mubarsq); }
            break;
        case 63: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},34,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},49,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},25,mubarsq); }
            else { return LFh63(m, mubarsq); }
            break;
        case 64: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},35,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},50,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},26,mubarsq); }
            else { return LFh64(m, mubarsq); }
            break;
        case 65: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},37,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},58,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},27,mubarsq); }
            else { return LFh65(m, mubarsq); }
            break;
        case 66: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},39,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},31,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},23,mubarsq); }
            else { return LFh66(m, mubarsq); }
            break;
        case 67: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},40,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},32,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},24,mubarsq); }
            else { return LFh67(m, mubarsq); }
            break;
        case 68: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},42,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},42,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},25,mubarsq); }
            else { return LFh68(m, mubarsq); }
            break;
        case 69: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},43,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},43,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},26,mubarsq); }
            else { return LFh69(m, mubarsq); }
            break;
        case 70: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},47,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},54,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},27,mubarsq); }
            else { return LFh70(m, mubarsq); }
            break;
        case 71: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},49,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},34,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},25,mubarsq); }
            else { return LFh71(m, mubarsq); }
            break;
        case 72: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},50,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},35,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},26,mubarsq); }
            else { return LFh72(m, mubarsq); }
            break;
        case 73: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},54,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},47,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},27,mubarsq); }
            else { return LFh73(m, mubarsq); }
            break;
        case 74: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},58,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},37,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},27,mubarsq); }
            else { return LFh74(m, mubarsq); }
            break;
        case 75: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},38,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},38,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},30,mubarsq); }
            else { return LFh75(m, mubarsq); }
            break;
        case 76: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},39,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},39,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},31,mubarsq); }
            else { return LFh76(m, mubarsq); }
            break;
        case 77: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},40,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},40,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},32,mubarsq); }
            else { return LFh77(m, mubarsq); }
            break;
        case 78: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},41,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},48,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},33,mubarsq); }
            else { return LFh78(m, mubarsq); }
            break;
        case 79: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},42,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},49,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},34,mubarsq); }
            else { return LFh79(m, mubarsq); }
            break;
        case 80: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},43,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},50,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},35,mubarsq); }
            else { return LFh80(m, mubarsq); }
            break;
        case 81: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},47,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},58,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},37,mubarsq); }
            else { return LFh81(m, mubarsq); }
            break;
        case 82: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},48,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},41,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},33,mubarsq); }
            else { return LFh82(m, mubarsq); }
            break;
        case 83: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},49,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},42,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},34,mubarsq); }
            else { return LFh83(m, mubarsq); }
            break;
        case 84: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},50,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},43,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},35,mubarsq); }
            else { return LFh84(m, mubarsq); }
            break;
        case 85: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},53,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},53,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},36,mubarsq); }
            else { return LFh85(m, mubarsq); }
            break;
        case 86: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},54,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},54,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},37,mubarsq); }
            else { return LFh86(m, mubarsq); }
            break;
        case 87: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},57,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},46,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},36,mubarsq); }
            else { return LFh87(m, mubarsq); }
            break;
        case 88: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},58,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},47,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},37,mubarsq); }
            else { return LFh88(m, mubarsq); }
            break;
        case 89: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},48,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},48,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},41,mubarsq); }
            else { return LFh89(m, mubarsq); }
            break;
        case 90: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},49,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},49,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},42,mubarsq); }
            else { return LFh90(m, mubarsq); }
            break;
        case 91: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},50,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},50,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},43,mubarsq); }
            else { return LFh91(m, mubarsq); }
            break;
        case 92: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},53,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},57,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},46,mubarsq); }
            else { return LFh92(m, mubarsq); }
            break;
        case 93: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},54,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},58,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},47,mubarsq); }
            else { return LFh93(m, mubarsq); }
            break;
        case 94: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},55,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},51,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},44,mubarsq); }
            else { return LFh94(m, mubarsq); }
            break;
        case 95: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},56,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},52,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},45,mubarsq); }
            else { return LFh95(m, mubarsq); }
            break;
        case 96: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},57,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},53,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},46,mubarsq); }
            else { return LFh96(m, mubarsq); }
            break;
        case 97: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},58,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},54,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},47,mubarsq); }
            else { return LFh97(m, mubarsq); }
            break;
        case 98: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},55,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},55,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},51,mubarsq); }
            else { return LFh98(m, mubarsq); }
            break;
        case 99: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},56,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},56,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},52,mubarsq); }
            else { return LFh99(m, mubarsq); }
            break;
        case 100: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},57,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},57,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},53,mubarsq); }
            else { return LFh100(m, mubarsq); }
            break;
        case 101: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2]},58,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1]},58,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1]},54,mubarsq); }
            else { return LFh101(m, mubarsq); }
            break;
        case 102: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},76,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},76,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},76,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},66,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},66,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},61,mubarsq); }
            else { return LFh102(m, mubarsq); }
            break;
        case 103: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},77,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},77,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},77,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},67,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},67,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},62,mubarsq); }
            else { return LFh103(m, mubarsq); }
            break;
        case 104: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},79,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},79,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},90,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},68,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},71,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},63,mubarsq); }
            else { return LFh104(m, mubarsq); }
            break;
        case 105: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},80,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},80,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},91,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},69,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},72,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},64,mubarsq); }
            else { return LFh105(m, mubarsq); }
            break;
        case 106: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},81,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},81,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},101,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},70,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},74,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},65,mubarsq); }
            else { return LFh106(m, mubarsq); }
            break;
        case 107: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},83,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},90,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},79,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},71,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},68,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},63,mubarsq); }
            else { return LFh107(m, mubarsq); }
            break;
        case 108: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},84,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},91,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},80,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},72,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},69,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},64,mubarsq); }
            else { return LFh108(m, mubarsq); }
            break;
        case 109: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},86,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},93,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},93,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},73,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},73,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},65,mubarsq); }
            else { return LFh109(m, mubarsq); }
            break;
        case 110: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},88,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},101,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},81,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},74,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},70,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},65,mubarsq); }
            else { return LFh110(m, mubarsq); }
            break;
        case 111: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},90,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},83,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},83,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},71,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},71,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},68,mubarsq); }
            else { return LFh111(m, mubarsq); }
            break;
        case 112: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},91,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},84,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},84,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},72,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},72,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},69,mubarsq); }
            else { return LFh112(m, mubarsq); }
            break;
        case 113: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},93,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},86,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},97,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},73,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},74,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},70,mubarsq); }
            else { return LFh113(m, mubarsq); }
            break;
        case 114: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},97,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},97,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},86,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},74,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},73,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},70,mubarsq); }
            else { return LFh114(m, mubarsq); }
            break;
        case 115: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},101,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},88,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},88,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},74,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},74,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},73,mubarsq); }
            else { return LFh115(m, mubarsq); }
            break;
        case 116: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},89,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},89,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},89,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},82,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},82,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},78,mubarsq); }
            else { return LFh116(m, mubarsq); }
            break;
        case 117: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},90,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},90,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},90,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},83,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},83,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},79,mubarsq); }
            else { return LFh117(m, mubarsq); }
            break;
        case 118: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},91,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},91,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},91,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},84,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},84,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},80,mubarsq); }
            else { return LFh118(m, mubarsq); }
            break;
        case 119: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},93,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},93,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},101,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},86,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},88,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},81,mubarsq); }
            else { return LFh119(m, mubarsq); }
            break;
        case 120: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},97,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},101,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},93,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},88,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},86,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},81,mubarsq); }
            else { return LFh120(m, mubarsq); }
            break;
        case 121: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},100,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},96,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},96,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},87,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},87,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},85,mubarsq); }
            else { return LFh121(m, mubarsq); }
            break;
        case 122: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},101,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},97,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},97,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},88,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},88,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},86,mubarsq); }
            else { return LFh122(m, mubarsq); }
            break;
        case 123: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},100,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},100,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},100,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},96,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},96,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},92,mubarsq); }
            else { return LFh123(m, mubarsq); }
            break;
        case 124: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3]},101,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},101,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},101,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3]},97,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},97,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2]},93,mubarsq); }
            else { return LFh124(m, mubarsq); }
            break;
        case 125: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3],m[4]},117,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},117,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},117,mubarsq); }
            else if (rel_diff(m[0],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},117,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},111,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},111,mubarsq); }
            else if (rel_diff(m[1],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},111,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},107,mubarsq); }
            else if (rel_diff(m[2],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},107,mubarsq); }
            else if (rel_diff(m[3],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},104,mubarsq); }
            else { return LFh125(m, mubarsq); }
            break;
        case 126: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3],m[4]},118,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},118,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},118,mubarsq); }
            else if (rel_diff(m[0],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},118,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},112,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},112,mubarsq); }
            else if (rel_diff(m[1],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},112,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},108,mubarsq); }
            else if (rel_diff(m[2],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},108,mubarsq); }
            else if (rel_diff(m[3],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},105,mubarsq); }
            else { return LFh126(m, mubarsq); }
            break;
        case 127: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3],m[4]},119,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},119,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},119,mubarsq); }
            else if (rel_diff(m[0],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},124,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},113,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},113,mubarsq); }
            else if (rel_diff(m[1],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},115,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},109,mubarsq); }
            else if (rel_diff(m[2],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},110,mubarsq); }
            else if (rel_diff(m[3],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},106,mubarsq); }
            else { return LFh127(m, mubarsq); }
            break;
        case 128: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3],m[4]},120,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},120,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},124,mubarsq); }
            else if (rel_diff(m[0],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},119,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},114,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},115,mubarsq); }
            else if (rel_diff(m[1],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},113,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},110,mubarsq); }
            else if (rel_diff(m[2],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},109,mubarsq); }
            else if (rel_diff(m[3],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},106,mubarsq); }
            else { return LFh128(m, mubarsq); }
            break;
        case 129: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3],m[4]},122,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},124,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},120,mubarsq); }
            else if (rel_diff(m[0],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},120,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},115,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},114,mubarsq); }
            else if (rel_diff(m[1],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},114,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},110,mubarsq); }
            else if (rel_diff(m[2],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},110,mubarsq); }
            else if (rel_diff(m[3],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},109,mubarsq); }
            else { return LFh129(m, mubarsq); }
            break;
        case 130: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3],m[4]},124,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},122,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},122,mubarsq); }
            else if (rel_diff(m[0],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},122,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},115,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},115,mubarsq); }
            else if (rel_diff(m[1],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},115,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},114,mubarsq); }
            else if (rel_diff(m[2],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},114,mubarsq); }
            else if (rel_diff(m[3],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},113,mubarsq); }
            else { return LFh130(m, mubarsq); }
            break;
        case 131: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3],m[4]},124,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},124,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},124,mubarsq); }
            else if (rel_diff(m[0],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},124,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4]},122,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},122,mubarsq); }
            else if (rel_diff(m[1],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},122,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4]},120,mubarsq); }
            else if (rel_diff(m[2],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},120,mubarsq); }
            else if (rel_diff(m[3],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3]},119,mubarsq); }
            else { return LFh131(m, mubarsq); }
            break;
        case 132: 
            if (rel_diff(m[0],m[1]) <= 5e-3) { return LF({m[0],m[2],m[3],m[4],m[5]},131,mubarsq); }
            else if (rel_diff(m[0],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4],m[5]},131,mubarsq); }
            else if (rel_diff(m[0],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4],m[5]},131,mubarsq); }
            else if (rel_diff(m[0],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[5]},131,mubarsq); }
            else if (rel_diff(m[0],m[5]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[4]},131,mubarsq); }
            else if (rel_diff(m[1],m[2]) <= 5e-3) { return LF({m[0],m[1],m[3],m[4],m[5]},130,mubarsq); }
            else if (rel_diff(m[1],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4],m[5]},130,mubarsq); }
            else if (rel_diff(m[1],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[5]},130,mubarsq); }
            else if (rel_diff(m[1],m[5]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[4]},130,mubarsq); }
            else if (rel_diff(m[2],m[3]) <= 5e-3) { return LF({m[0],m[1],m[2],m[4],m[5]},129,mubarsq); }
            else if (rel_diff(m[2],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[5]},129,mubarsq); }
            else if (rel_diff(m[2],m[5]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[4]},129,mubarsq); }
            else if (rel_diff(m[3],m[4]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[5]},128,mubarsq); }
            else if (rel_diff(m[3],m[5]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[4]},128,mubarsq); }
            else if (rel_diff(m[4],m[5]) <= 5e-3) { return LF({m[0],m[1],m[2],m[3],m[4]},127,mubarsq); }
            else { return LFh132(m, mubarsq); }
            break;
    }
    return 0.0;
}
