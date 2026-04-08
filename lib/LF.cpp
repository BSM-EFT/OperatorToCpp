#include "pch.h"
#include "LF.h"
#include "complex_math.h"
using std::vector;
using std::complex;

// elaborate definitions of loop-functions in terms of masses
complex<double> LF(vector<complex<double> > masses, int code, double mubarsq) {
    switch(code) {
        case 1: 
            return pow(masses[0],2)*(1 + log(mubarsq/pow(masses[0],2)));
            break;
        case 2: 
            return 1 + log(mubarsq/pow(masses[0],2));
            break;
        case 3: 
            return (1 + log(mubarsq/pow(masses[0],2)))/pow(masses[0],2);
            break;
        case 4: 
            return pow(masses[0],2)*(1 + 2*log(mubarsq/pow(masses[0],2)));
            break;
        case 5: 
            return log(mubarsq/pow(masses[0],2));
            break;
        case 6: 
            return -pow(masses[0],-2);
            break;
        case 7: 
            return -0.5 + log(mubarsq/pow(masses[0],2));
            break;
        case 8: 
            return -0.5*1/pow(masses[0],2);
            break;
        case 9: 
            return -0.8333333333333334 + log(mubarsq/pow(masses[0],2));
            break;
        case 10: 
            return -0.3333333333333333*1/pow(masses[0],2);
            break;
        case 11: 
            return 1/(6*pow(masses[0],4));
            break;
        case 12: 
            return -0.25*1/pow(masses[0],2);
            break;
        case 13: 
            return 1/(12*pow(masses[0],4));
            break;
        case 14: 
            return -0.08333333333333333*1/pow(masses[0],6);
            break;
        case 15: 
            return -0.2*1/pow(masses[0],2);
            break;
        case 16: 
            return 1/(20*pow(masses[0],4));
            break;
        case 17: 
            return -0.03333333333333333*1/pow(masses[0],6);
            break;
        case 18: 
            return 1/(20*pow(masses[0],8));
            break;
        case 19: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},4,mubarsq); }
            else { return (pow(masses[0],4) - pow(masses[1],4) + pow(masses[0],4)*log(mubarsq/pow(masses[0],2)) - pow(masses[1],4)*log(mubarsq/pow(masses[1],2)))/(pow(masses[0],2) - pow(masses[1],2)); }
            break;
        case 20: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},5,mubarsq); }
            else { return (pow(masses[0],2) - pow(masses[1],2) + pow(masses[0],2)*log(mubarsq/pow(masses[0],2)) - pow(masses[1],2)*log(mubarsq/pow(masses[1],2)))/(pow(masses[0],2) - pow(masses[1],2)); }
            break;
        case 21: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},7,mubarsq); }
            else { return (pow(masses[0],4) - pow(masses[0],2)*pow(masses[1],2) + pow(masses[0],4)*log(mubarsq/pow(masses[0],2)) + (-2*pow(masses[0],2)*pow(masses[1],2) + pow(masses[1],4))*log(mubarsq/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),2); }
            break;
        case 22: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},8,mubarsq); }
            else { return (pow(masses[0],2) - pow(masses[1],2) - pow(masses[0],2)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),2); }
            break;
        case 23: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},10,mubarsq); }
            else { return (3*pow(masses[0],4) - 4*pow(masses[0],2)*pow(masses[1],2) + pow(masses[1],4) - 2*pow(masses[0],4)*log(pow(masses[0],2)/pow(masses[1],2)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),3)); }
            break;
        case 24: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},11,mubarsq); }
            else { return (-pow(masses[0],4) + pow(masses[1],4) + 2*pow(masses[0],2)*pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/(2*pow(masses[1],2)*pow(-pow(masses[0],2) + pow(masses[1],2),3)); }
            break;
        case 25: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},13,mubarsq); }
            else { return (2*pow(masses[0],6) + 3*pow(masses[0],4)*pow(masses[1],2) - 6*pow(masses[0],2)*pow(masses[1],4) + pow(masses[1],6) - 6*pow(masses[0],4)*pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/(6*pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[1],2),4)); }
            break;
        case 26: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},14,mubarsq); }
            else { return -0.16666666666666666*(pow(masses[0],6) - 6*pow(masses[0],4)*pow(masses[1],2) + 3*pow(masses[0],2)*pow(masses[1],4) + 2*pow(masses[1],6) + 6*pow(masses[0],2)*pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(-(pow(masses[0],2)*masses[1]) + pow(masses[1],3),4); }
            break;
        case 27: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},18,mubarsq); }
            else { return (-pow(masses[0],8) + 6*pow(masses[0],6)*pow(masses[1],2) - 18*pow(masses[0],4)*pow(masses[1],4) + 10*pow(masses[0],2)*pow(masses[1],6) + 3*pow(masses[1],8) + 12*pow(masses[0],2)*pow(masses[1],6)*log(pow(masses[0],2)/pow(masses[1],2)))/(12*pow(masses[1],6)*pow(-pow(masses[0],2) + pow(masses[1],2),5)); }
            break;
        case 28: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},7,mubarsq); }
            else { return (-(pow(masses[0],2)*pow(masses[1],2)) + pow(masses[1],4) + (pow(masses[0],4) - 2*pow(masses[0],2)*pow(masses[1],2))*log(mubarsq/pow(masses[0],2)) + pow(masses[1],4)*log(mubarsq/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),2); }
            break;
        case 29: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},8,mubarsq); }
            else { return (-pow(masses[0],2) + pow(masses[1],2) + pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),2); }
            break;
        case 30: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},9,mubarsq); }
            else { return (-2*pow(masses[0],4)*pow(masses[1],2) + 2*pow(masses[0],2)*pow(masses[1],4) + (pow(masses[0],6) - 3*pow(masses[0],4)*pow(masses[1],2))*log(mubarsq/pow(masses[0],2)) + (3*pow(masses[0],2)*pow(masses[1],4) - pow(masses[1],6))*log(mubarsq/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),3); }
            break;
        case 31: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},10,mubarsq); }
            else { return (-pow(masses[0],4) + pow(masses[1],4) + 2*pow(masses[0],2)*pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),3); }
            break;
        case 32: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},11,mubarsq); }
            else { return (-2*pow(masses[0],2) + 2*pow(masses[1],2) + (pow(masses[0],2) + pow(masses[1],2))*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),3); }
            break;
        case 33: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},12,mubarsq); }
            else { return -0.5*(2*pow(masses[0],6) + 3*pow(masses[0],4)*pow(masses[1],2) - 6*pow(masses[0],2)*pow(masses[1],4) + pow(masses[1],6) - 6*pow(masses[0],4)*pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),4); }
            break;
        case 34: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},13,mubarsq); }
            else { return (-5*pow(masses[0],4) + 4*pow(masses[0],2)*pow(masses[1],2) + pow(masses[1],4) + 2*(pow(masses[0],4) + 2*pow(masses[0],2)*pow(masses[1],2))*log(pow(masses[0],2)/pow(masses[1],2)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),4)); }
            break;
        case 35: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},14,mubarsq); }
            else { return -0.5*(pow(masses[0],4) + 4*pow(masses[0],2)*pow(masses[1],2) - 5*pow(masses[1],4) - 2*(2*pow(masses[0],2)*pow(masses[1],2) + pow(masses[1],4))*log(pow(masses[0],2)/pow(masses[1],2)))/(pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[1],2),4)); }
            break;
        case 36: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},17,mubarsq); }
            else { return (pow(masses[0],6) + 9*pow(masses[0],4)*pow(masses[1],2) - 9*pow(masses[0],2)*pow(masses[1],4) - pow(masses[1],6) - 6*pow(masses[0],2)*pow(masses[1],2)*(pow(masses[0],2) + pow(masses[1],2))*log(pow(masses[0],2)/pow(masses[1],2)))/(3*pow(masses[1],2)*pow(-pow(masses[0],2) + pow(masses[1],2),5)); }
            break;
        case 37: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},18,mubarsq); }
            else { return -0.16666666666666666*(pow(masses[0],6) - 9*pow(masses[0],4)*pow(masses[1],2) - 9*pow(masses[0],2)*pow(masses[1],4) + 17*pow(masses[1],6) + 6*(3*pow(masses[0],2)*pow(masses[1],4) + pow(masses[1],6))*log(pow(masses[0],2)/pow(masses[1],2)))/(pow(masses[1],4)*pow(-pow(masses[0],2) + pow(masses[1],2),5)); }
            break;
        case 38: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},9,mubarsq); }
            else { return (-pow(masses[0],6) + 2*pow(masses[0],4)*pow(masses[1],2) + pow(masses[0],2)*pow(masses[1],4) - 2*pow(masses[1],6) + 2*(pow(masses[0],6) - 3*pow(masses[0],4)*pow(masses[1],2) + 3*pow(masses[0],2)*pow(masses[1],4))*log(mubarsq/pow(masses[0],2)) - 2*pow(masses[1],6)*log(mubarsq/pow(masses[1],2)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),3)); }
            break;
        case 39: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},10,mubarsq); }
            else { return -0.5*(pow(masses[0],4) - 4*pow(masses[0],2)*pow(masses[1],2) + 3*pow(masses[1],4) + 2*pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),3); }
            break;
        case 40: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},11,mubarsq); }
            else { return (pow(masses[0],4) - pow(masses[1],4) - 2*pow(masses[0],2)*pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/(2*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),3)); }
            break;
        case 41: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},12,mubarsq); }
            else { return -0.5*(pow(masses[0],6) - 6*pow(masses[0],4)*pow(masses[1],2) + 3*pow(masses[0],2)*pow(masses[1],4) + 2*pow(masses[1],6) + 6*pow(masses[0],2)*pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),4); }
            break;
        case 42: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},13,mubarsq); }
            else { return (pow(masses[0],4) + 4*pow(masses[0],2)*pow(masses[1],2) - 5*pow(masses[1],4) - 2*(2*pow(masses[0],2)*pow(masses[1],2) + pow(masses[1],4))*log(pow(masses[0],2)/pow(masses[1],2)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),4)); }
            break;
        case 43: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},14,mubarsq); }
            else { return -0.5*(-5*pow(masses[0],4) + 4*pow(masses[0],2)*pow(masses[1],2) + pow(masses[1],4) + 2*(pow(masses[0],4) + 2*pow(masses[0],2)*pow(masses[1],2))*log(pow(masses[0],2)/pow(masses[1],2)))/(pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),4)); }
            break;
        case 44: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},15,mubarsq); }
            else { return (-pow(masses[0],8) + 8*pow(masses[0],6)*pow(masses[1],2) - 8*pow(masses[0],2)*pow(masses[1],6) + pow(masses[1],8) - 12*pow(masses[0],4)*pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 45: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},16,mubarsq); }
            else { return (pow(masses[0],6) + 9*pow(masses[0],4)*pow(masses[1],2) - 9*pow(masses[0],2)*pow(masses[1],4) - pow(masses[1],6) - 6*pow(masses[0],2)*pow(masses[1],2)*(pow(masses[0],2) + pow(masses[1],2))*log(pow(masses[0],2)/pow(masses[1],2)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 46: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},17,mubarsq); }
            else { return (3*(pow(masses[0],4) - pow(masses[1],4)) - (pow(masses[0],4) + 4*pow(masses[0],2)*pow(masses[1],2) + pow(masses[1],4))*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),5); }
            break;
        case 47: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},18,mubarsq); }
            else { return (pow(masses[0],6) + 9*pow(masses[0],4)*pow(masses[1],2) - 9*pow(masses[0],2)*pow(masses[1],4) - pow(masses[1],6) - 6*pow(masses[0],2)*pow(masses[1],2)*(pow(masses[0],2) + pow(masses[1],2))*log(pow(masses[0],2)/pow(masses[1],2)))/(2*pow(masses[0],2)*pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 48: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},12,mubarsq); }
            else { return (-2*pow(masses[0],6) + 9*pow(masses[0],4)*pow(masses[1],2) - 18*pow(masses[0],2)*pow(masses[1],4) + 11*pow(masses[1],6) + 6*pow(masses[1],6)*log(pow(masses[0],2)/pow(masses[1],2)))/(6*pow(pow(masses[0],2) - pow(masses[1],2),4)); }
            break;
        case 49: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},13,mubarsq); }
            else { return (pow(masses[0],6) - 6*pow(masses[0],4)*pow(masses[1],2) + 3*pow(masses[0],2)*pow(masses[1],4) + 2*pow(masses[1],6) + 6*pow(masses[0],2)*pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/(6*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),4)); }
            break;
        case 50: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},14,mubarsq); }
            else { return -0.16666666666666666*(2*pow(masses[0],6) + 3*pow(masses[0],4)*pow(masses[1],2) - 6*pow(masses[0],2)*pow(masses[1],4) + pow(masses[1],6) - 6*pow(masses[0],4)*pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],3) - masses[0]*pow(masses[1],2),4); }
            break;
        case 51: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},15,mubarsq); }
            else { return (-pow(masses[0],8) + 6*pow(masses[0],6)*pow(masses[1],2) - 18*pow(masses[0],4)*pow(masses[1],4) + 10*pow(masses[0],2)*pow(masses[1],6) + 3*pow(masses[1],8) + 12*pow(masses[0],2)*pow(masses[1],6)*log(pow(masses[0],2)/pow(masses[1],2)))/(3*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 52: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},16,mubarsq); }
            else { return (pow(masses[0],6) - 9*pow(masses[0],4)*pow(masses[1],2) - 9*pow(masses[0],2)*pow(masses[1],4) + 17*pow(masses[1],6) + 6*(3*pow(masses[0],2)*pow(masses[1],4) + pow(masses[1],6))*log(pow(masses[0],2)/pow(masses[1],2)))/(6*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 53: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},17,mubarsq); }
            else { return (-pow(masses[0],6) - 9*pow(masses[0],4)*pow(masses[1],2) + 9*pow(masses[0],2)*pow(masses[1],4) + pow(masses[1],6) + 6*pow(masses[0],2)*pow(masses[1],2)*(pow(masses[0],2) + pow(masses[1],2))*log(pow(masses[0],2)/pow(masses[1],2)))/(3*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 54: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},18,mubarsq); }
            else { return (-17*pow(masses[0],6) + 9*pow(masses[0],4)*pow(masses[1],2) + 9*pow(masses[0],2)*pow(masses[1],4) - pow(masses[1],6) + 6*(pow(masses[0],6) + 3*pow(masses[0],4)*pow(masses[1],2))*log(pow(masses[0],2)/pow(masses[1],2)))/(6*pow(masses[0],4)*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 55: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},15,mubarsq); }
            else { return -0.08333333333333333*(3*pow(masses[0],8) - 16*pow(masses[0],6)*pow(masses[1],2) + 36*pow(masses[0],4)*pow(masses[1],4) - 48*pow(masses[0],2)*pow(masses[1],6) + 25*pow(masses[1],8) + 12*pow(masses[1],8)*log(pow(masses[0],2)/pow(masses[1],2)))/pow(pow(masses[0],2) - pow(masses[1],2),5); }
            break;
        case 56: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},16,mubarsq); }
            else { return (pow(masses[0],8) - 6*pow(masses[0],6)*pow(masses[1],2) + 18*pow(masses[0],4)*pow(masses[1],4) - 10*pow(masses[0],2)*pow(masses[1],6) - 3*pow(masses[1],8) - 12*pow(masses[0],2)*pow(masses[1],6)*log(pow(masses[0],2)/pow(masses[1],2)))/(12*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 57: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},17,mubarsq); }
            else { return (-pow(masses[0],8) + 8*pow(masses[0],6)*pow(masses[1],2) - 8*pow(masses[0],2)*pow(masses[1],6) + pow(masses[1],8) - 12*pow(masses[0],4)*pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/(12*pow(masses[0],4)*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 58: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0]},18,mubarsq); }
            else { return -0.08333333333333333*(-3*pow(masses[0],8) - 10*pow(masses[0],6)*pow(masses[1],2) + 18*pow(masses[0],4)*pow(masses[1],4) - 6*pow(masses[0],2)*pow(masses[1],6) + pow(masses[1],8) + 12*pow(masses[0],6)*pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/(pow(masses[0],6)*pow(pow(masses[0],2) - pow(masses[1],2),5)); }
            break;
        case 59: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},28,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},28,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},21,mubarsq); }
            else { return (pow(masses[0],4)*(pow(masses[1],2) - pow(masses[2],2))*log(mubarsq/pow(masses[0],2)) + pow(masses[1],4)*(-pow(masses[0],2) + pow(masses[2],2))*log(mubarsq/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2)) + pow(masses[2],4)*log(mubarsq/pow(masses[2],2))))/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 60: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},29,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},29,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},22,mubarsq); }
            else { return (pow(masses[1],2)*(-pow(masses[0],2) + pow(masses[2],2))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2)))/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 61: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},31,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},39,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},23,mubarsq); }
            else { return ((pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/(-pow(masses[0],2) + pow(masses[1],2)) + (pow(masses[2],2)*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2)) - (pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(-2*pow(masses[1],2) + pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[2],2))))/pow(pow(masses[0],2) - pow(masses[2],2),2))/pow(pow(masses[1],2) - pow(masses[2],2),2); }
            break;
        case 62: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},32,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},40,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},24,mubarsq); }
            else { return ((pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/(-pow(masses[0],2) + pow(masses[1],2)) + ((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2)) + (pow(masses[0],2)*pow(masses[1],2) - pow(masses[2],4))*log(pow(masses[0],2)/pow(masses[2],2)))/pow(pow(masses[0],2) - pow(masses[2],2),2))/pow(pow(masses[1],2) - pow(masses[2],2),2); }
            break;
        case 63: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},34,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},49,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},25,mubarsq); }
            else { return ((2*pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/(-pow(masses[0],2) + pow(masses[1],2)) + (-((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[0],2)*(-3*pow(masses[1],2) + pow(masses[2],2)) + pow(masses[2],2)*(pow(masses[1],2) + pow(masses[2],2)))) + 2*(pow(masses[0],4)*pow(masses[1],4) + pow(masses[1],2)*pow(masses[2],6) + pow(masses[0],2)*(-3*pow(masses[1],2)*pow(masses[2],4) + pow(masses[2],6)))*log(pow(masses[0],2)/pow(masses[2],2)))/pow(pow(masses[0],2) - pow(masses[2],2),3))/(2*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 64: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},35,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},50,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},26,mubarsq); }
            else { return (2*pow(masses[1],2)*pow(masses[2],2)*pow(-pow(masses[0],2) + pow(masses[2],2),3)*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[2],2)*(pow(masses[1],2) - 3*pow(masses[2],2)) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))) - 2*pow(masses[2],2)*(pow(masses[0],4)*pow(masses[1],2) + pow(masses[2],6) + pow(masses[0],2)*(pow(masses[1],4) - 3*pow(masses[1],2)*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[2],2))))/(2*(pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 65: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},37,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},58,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},27,mubarsq); }
            else { return (-6*pow(masses[1],2)*pow(-(pow(masses[0],2)*masses[2]) + pow(masses[2],3),4)*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(2*pow(masses[1],4)*pow(masses[2],4) - 7*pow(masses[1],2)*pow(masses[2],6) + 11*pow(masses[2],8) + pow(masses[0],4)*(-pow(masses[1],4) + 5*pow(masses[1],2)*pow(masses[2],2) + 2*pow(masses[2],4)) + pow(masses[0],2)*(5*pow(masses[1],4)*pow(masses[2],2) - 10*pow(masses[1],2)*pow(masses[2],4) - 7*pow(masses[2],6))) + 6*pow(masses[2],4)*(pow(masses[0],6)*pow(masses[1],2) - pow(masses[2],8) + pow(masses[0],4)*(pow(masses[1],4) - 4*pow(masses[1],2)*pow(masses[2],2)) + pow(masses[0],2)*(pow(masses[1],6) - 4*pow(masses[1],4)*pow(masses[2],2) + 6*pow(masses[1],2)*pow(masses[2],4)))*log(pow(masses[0],2)/pow(masses[2],2))))/(6*(pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],4)*pow(pow(masses[0],2) - pow(masses[2],2),4)*pow(pow(masses[1],2) - pow(masses[2],2),4)); }
            break;
        case 66: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},39,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},31,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},23,mubarsq); }
            else { return (pow(masses[1],2)*(-2*pow(masses[0],2)*pow(masses[2],4) + pow(masses[1],2)*pow(masses[2],4) - pow(masses[0],4)*(pow(masses[1],2) - 2*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*(pow(masses[1],2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2)) + (-pow(masses[0],2) + pow(masses[1],2))*pow(masses[2],4)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 67: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},40,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},32,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},24,mubarsq); }
            else { return ((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],4) + pow(masses[0],2)*pow(masses[2],2))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2)) + (-pow(masses[0],2) + pow(masses[1],2))*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 68: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},42,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},42,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},25,mubarsq); }
            else { return (-(pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[1],4) - 2*pow(masses[0],2)*pow(masses[2],2) + pow(masses[1],2)*pow(masses[2],2))*log(pow(masses[0],2)/pow(masses[1],2))) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(-2*pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))) - (pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],2)*(2*pow(masses[0],2)*pow(masses[1],2) - pow(masses[2],2)*(pow(masses[1],2) + pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 69: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},43,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},43,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},26,mubarsq); }
            else { return (pow(pow(masses[0],2) - pow(masses[2],2),2)*(-2*pow(masses[1],4) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(2*pow(masses[0],2) - pow(masses[1],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2)) + (pow(masses[0],2) - pow(masses[1],2))*(-2*pow(masses[2],4) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 70: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},47,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},54,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},27,mubarsq); }
            else { return pow(masses[0],2)/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),3)) + 1/(2*(pow(masses[0],2) - pow(masses[2],2))*pow(-(pow(masses[1],2)*masses[2]) + pow(masses[2],3),2)) + log(pow(masses[0],2)/pow(masses[1],2))/((-pow(masses[0],2) + pow(masses[1],2))*pow(pow(masses[1],2) - pow(masses[2],2),3)) - (pow(masses[1],2)*(-3*pow(masses[0],2) + 4*pow(masses[1],2) - pow(masses[2],2))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),4)) + ((2*pow(masses[0],2) + pow(masses[1],2) - 3*pow(masses[2],2))*log(pow(masses[0],2)/pow(masses[2],2)))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(-pow(masses[1],2) + pow(masses[2],2),3)) - (pow(masses[2],2)*(3*pow(masses[0],4) + pow(masses[1],4) - 4*pow(masses[1],2)*pow(masses[2],2) + 6*pow(masses[2],4) + 2*pow(masses[0],2)*(pow(masses[1],2) - 4*pow(masses[2],2)))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),4)); }
            break;
        case 71: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},49,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},34,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},25,mubarsq); }
            else { return (-2*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],6)*pow(masses[2],2) + pow(masses[0],4)*pow(masses[2],4) + pow(masses[0],2)*(pow(masses[1],6) - 3*pow(masses[1],4)*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[0],2)*(pow(masses[1],2) - 3*pow(masses[2],2)) + pow(masses[1],2)*(pow(masses[1],2) + pow(masses[2],2))) + 2*pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],4)*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(pow(masses[0],2) - pow(masses[1],2),3)*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 72: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},50,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},35,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},26,mubarsq); }
            else { return (2*pow(masses[1],2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],6) - 3*pow(masses[0],2)*pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*pow(masses[2],2)*(pow(masses[0],2) + pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(-3*pow(masses[1],4) + pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))) - 2*pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(masses[1],2)*pow(-pow(masses[0],2) + pow(masses[1],2),3)*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 73: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},54,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},47,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},27,mubarsq); }
            else { return pow(masses[0],2)/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),2)) + 1/(2*(pow(masses[0],2) - pow(masses[1],2))*pow(pow(masses[1],3) - masses[1]*pow(masses[2],2),2)) + ((2*pow(masses[0],2) - 3*pow(masses[1],2) + pow(masses[2],2))*log(pow(masses[0],2)/pow(masses[1],2)))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),3)) + (pow(masses[1],2)*(3*pow(masses[0],4) + 6*pow(masses[1],4) - 4*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4) + pow(masses[0],2)*(-8*pow(masses[1],2) + 2*pow(masses[2],2)))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(-pow(masses[0],2) + pow(masses[1],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),4)) + log(pow(masses[0],2)/pow(masses[2],2))/((-pow(masses[0],2) + pow(masses[2],2))*pow(-pow(masses[1],2) + pow(masses[2],2),3)) - (pow(masses[2],2)*(-3*pow(masses[0],2) - pow(masses[1],2) + 4*pow(masses[2],2))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),4)); }
            break;
        case 74: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},58,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},37,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},27,mubarsq); }
            else { return (6*pow(masses[1],4)*(pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],8) + 6*pow(masses[0],2)*pow(masses[1],4)*pow(masses[2],2) - 4*pow(masses[0],2)*pow(masses[1],2)*pow(masses[2],2)*(pow(masses[0],2) + pow(masses[2],2)) + pow(masses[0],2)*pow(masses[2],2)*(pow(masses[0],4) + pow(masses[0],2)*pow(masses[2],2) + pow(masses[2],4)))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(-11*pow(masses[1],8) + 7*pow(masses[1],6)*pow(masses[2],2) - 2*pow(masses[1],4)*pow(masses[2],4) + pow(masses[0],4)*(-2*pow(masses[1],4) - 5*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4)) + pow(masses[0],2)*(7*pow(masses[1],6) + 10*pow(masses[1],4)*pow(masses[2],2) - 5*pow(masses[1],2)*pow(masses[2],4))) + 6*pow(masses[1],4)*pow(-pow(masses[0],2) + pow(masses[1],2),3)*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(6*pow(masses[1],4)*pow(pow(masses[0],2) - pow(masses[1],2),4)*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[2],2),4)); }
            break;
        case 75: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},38,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},38,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},30,mubarsq); }
            else { return (pow(masses[0],4)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[0],4) + 3*pow(masses[1],2)*pow(masses[2],2) - 2*pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2)))*log(mubarsq/pow(masses[0],2)) + pow(masses[1],6)*pow(pow(masses[0],2) - pow(masses[2],2),2)*log(mubarsq/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(-(pow(masses[1],2)*pow(masses[2],2)) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))) + (pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],6)*log(mubarsq/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 76: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},39,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},39,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},31,mubarsq); }
            else { return (pow(masses[1],4)*pow(pow(masses[0],2) - pow(masses[2],2),2)*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2)) + (pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],4)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 77: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},40,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},40,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},32,mubarsq); }
            else { return (pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2)) + (pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 78: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},41,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},48,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},33,mubarsq); }
            else { return (pow(masses[1],6)*pow(pow(masses[0],2) - pow(masses[2],2),3)*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[0],2)*pow(masses[2],4) - pow(masses[1],2)*pow(masses[2],4) + pow(masses[0],4)*(-pow(masses[1],2) + pow(masses[2],2))) + (pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],4)*(-(pow(masses[1],2)*pow(masses[2],2)) + pow(masses[0],2)*(3*pow(masses[1],2) - 2*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 79: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},42,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},49,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},34,mubarsq); }
            else { return (pow(masses[1],4)*pow(pow(masses[0],2) - pow(masses[2],2),3)*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(pow(masses[1],2) - 2*pow(masses[2],2))) + (pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],2)*(-pow(masses[2],4) + pow(masses[0],2)*(2*pow(masses[1],2) - pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 80: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},43,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},50,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},35,mubarsq); }
            else { return (pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - 2*pow(masses[1],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2)) + (pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2)*pow(masses[1],2) + pow(masses[1],2)*pow(masses[2],2) - 2*pow(masses[2],4))*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 81: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},47,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},58,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},37,mubarsq); }
            else { return (2*pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)) + ((pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(-5*pow(masses[1],4)*pow(masses[2],2) + 9*pow(masses[1],2)*pow(masses[2],4) - 2*pow(masses[2],6) + pow(masses[0],4)*(pow(masses[1],2) + pow(masses[2],2)) - pow(masses[0],2)*(pow(masses[1],4) - 2*pow(masses[1],2)*pow(masses[2],2) + 5*pow(masses[2],4))) - 2*(pow(masses[0],2) - pow(masses[1],2))*pow(masses[2],2)*(pow(masses[0],4)*pow(masses[1],2) + 2*pow(masses[0],2)*(pow(masses[1],4) - 2*pow(masses[1],2)*pow(masses[2],2)) + pow(masses[2],2)*(pow(masses[1],4) - 3*pow(masses[1],2)*pow(masses[2],2) + 3*pow(masses[2],4)))*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(masses[2],2)*pow(pow(masses[0],2) - pow(masses[2],2),4)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 82: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},48,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},41,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},33,mubarsq); }
            else { return (pow(masses[1],4)*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(2*pow(masses[1],2) - 3*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(-(pow(masses[0],2)*pow(masses[1],4)) + pow(masses[1],4)*pow(masses[2],2) + pow(masses[0],4)*(-pow(masses[1],2) + pow(masses[2],2))) - pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],6)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 83: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},49,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},42,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},34,mubarsq); }
            else { return (pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[1],4) + pow(masses[0],2)*(pow(masses[1],2) - 2*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(-2*pow(masses[1],2) + pow(masses[2],2))) - pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],4)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 84: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},50,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},43,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},35,mubarsq); }
            else { return -((pow(pow(masses[0],2) - pow(masses[2],2),2)*(-2*pow(masses[1],4) + pow(masses[0],2)*pow(masses[2],2) + pow(masses[1],2)*pow(masses[2],2))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) + pow(masses[1],2) - 2*pow(masses[2],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2)) - pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2))); }
            break;
        case 85: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},53,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},53,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},36,mubarsq); }
            else { return (2*pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*(pow(masses[1],4) - pow(masses[0],2)*pow(masses[2],2))*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[0],4)*(pow(masses[1],2) + pow(masses[2],2)) + pow(masses[1],2)*pow(masses[2],2)*(pow(masses[1],2) + pow(masses[2],2)) + pow(masses[0],2)*(pow(masses[1],4) - 6*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4))) - 2*pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],2)*(pow(masses[0],2)*pow(masses[1],2) - pow(masses[2],4))*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 86: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},54,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},54,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},37,mubarsq); }
            else { return (-(pow(pow(masses[0],2) - pow(masses[2],2),3)*(-3*pow(masses[1],4) + pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2))) + (pow(masses[0],2) - pow(masses[1],2))*(2*(pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[0],4) + pow(masses[1],4) - pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4) - pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))) + pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[2],2)*(pow(masses[1],2) - 3*pow(masses[2],2)) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 87: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},57,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},46,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},36,mubarsq); }
            else { return (2*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[1],8) + pow(masses[0],4)*pow(masses[2],4) + 2*pow(masses[0],2)*pow(masses[1],2)*(pow(masses[1],4) - 3*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4)))*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(-3*pow(masses[1],4)*pow(masses[2],2) + pow(masses[1],2)*pow(masses[2],4) + pow(masses[0],4)*(pow(masses[1],2) - 3*pow(masses[2],2)) + pow(masses[0],2)*(5*pow(masses[1],4) - 6*pow(masses[1],2)*pow(masses[2],2) + 5*pow(masses[2],4))) + 2*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],4)*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 88: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},58,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},47,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},37,mubarsq); }
            else { return (2*(3*pow(masses[1],6) - 3*pow(masses[1],4)*pow(masses[2],2) + pow(masses[0],2)*pow(masses[2],2)*(pow(masses[0],2) + 2*pow(masses[2],2)) + pow(masses[1],2)*(-4*pow(masses[0],2)*pow(masses[2],2) + pow(masses[2],4)))*log(pow(masses[0],2)/pow(masses[1],2)) + ((pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(-2*pow(masses[1],6) + 9*pow(masses[1],4)*pow(masses[2],2) - 5*pow(masses[1],2)*pow(masses[2],4) + pow(masses[0],4)*(pow(masses[1],2) + pow(masses[2],2)) - pow(masses[0],2)*(5*pow(masses[1],4) - 2*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4))) + 2*pow(masses[1],2)*pow(-pow(masses[0],2) + pow(masses[1],2),3)*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),2)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[1],2) - pow(masses[2],2),3)); }
            break;
        case 89: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},48,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},48,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},41,mubarsq); }
            else { return -0.5*(2*pow(masses[1],6)*pow(pow(masses[0],2) - pow(masses[2],2),3)*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[0],4) + 5*pow(masses[1],2)*pow(masses[2],2) - 3*pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))) - 2*pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],6)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),3)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 90: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},49,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},49,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},42,mubarsq); }
            else { return (-2*pow(masses[1],4)*pow(pow(masses[0],2) - pow(masses[2],2),3)*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[0],4) - 3*pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))) + 2*pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],4)*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),3)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 91: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},50,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},50,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},43,mubarsq); }
            else { return (-2*pow(masses[0],2)*pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(3*pow(masses[0],4) - pow(masses[1],2)*pow(masses[2],2) - pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))) + 2*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),3)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 92: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},53,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},57,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},46,mubarsq); }
            else { return (-2*pow(masses[1],4)*pow(pow(masses[0],2) - pow(masses[2],2),4)*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(5*pow(masses[1],4)*pow(masses[2],2) - 3*pow(masses[1],2)*pow(masses[2],4) + pow(masses[0],4)*(-3*pow(masses[1],2) + 5*pow(masses[2],2)) + pow(masses[0],2)*(pow(masses[1],4) - 6*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4))) + 2*pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(masses[2],2)*(pow(masses[2],2)*(pow(masses[1],2) - 2*pow(masses[2],2)) + pow(masses[0],2)*(2*pow(masses[1],2) - pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),4)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 93: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},54,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},58,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},47,mubarsq); }
            else { return (-2*pow(masses[0],2)*pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),4)*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(2*pow(masses[0],6) + pow(masses[1],2)*pow(masses[2],2)*(pow(masses[1],2) - pow(masses[2],2)) + pow(masses[0],4)*(-9*pow(masses[1],2) + 5*pow(masses[2],2)) + pow(masses[0],2)*(5*pow(masses[1],4) - 2*pow(masses[1],2)*pow(masses[2],2) - pow(masses[2],4))) + 2*pow(pow(masses[0],3) - masses[0]*pow(masses[1],2),2)*(pow(masses[0],2)*pow(masses[1],2) + 2*pow(masses[1],2)*pow(masses[2],2) - 3*pow(masses[2],4))*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),4)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 94: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},55,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},51,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},44,mubarsq); }
            else { return -0.5*(2*pow(masses[1],6)*pow(pow(masses[0],2) - pow(masses[2],2),3)*(pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(3*pow(masses[1],2) - 4*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(-4*pow(masses[0],2)*pow(masses[1],6)*pow(masses[2],2) + 2*pow(masses[1],6)*pow(masses[2],4) + pow(masses[0],8)*(-pow(masses[1],2) + pow(masses[2],2)) + pow(masses[0],6)*(5*pow(masses[1],4) - 2*pow(masses[1],2)*pow(masses[2],2) - 3*pow(masses[2],4)) + pow(masses[0],4)*(2*pow(masses[1],6) - 7*pow(masses[1],4)*pow(masses[2],2) + 7*pow(masses[1],2)*pow(masses[2],4))) + 2*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],8)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 95: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},56,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},52,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},45,mubarsq); }
            else { return (-2*pow(masses[1],4)*pow(pow(masses[0],2) - pow(masses[2],2),3)*(pow(masses[1],4) + pow(masses[0],2)*(2*pow(masses[1],2) - 3*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(-2*pow(masses[1],4)*pow(masses[2],4) + pow(masses[0],6)*(-pow(masses[1],2) + pow(masses[2],2)) + pow(masses[0],4)*(-5*pow(masses[1],4) + 2*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4)) + pow(masses[0],2)*(9*pow(masses[1],4)*pow(masses[2],2) - 5*pow(masses[1],2)*pow(masses[2],4))) - 2*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],6)*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 96: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},57,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},53,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},46,mubarsq); }
            else { return (-2*pow(masses[1],2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*(2*pow(masses[1],4) - pow(masses[1],2)*pow(masses[2],2) + pow(masses[0],2)*(pow(masses[1],2) - 2*pow(masses[2],2)))*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(-3*pow(masses[1],4)*pow(masses[2],2) + 5*pow(masses[1],2)*pow(masses[2],4) + pow(masses[0],4)*(5*pow(masses[1],2) - 3*pow(masses[2],2)) + pow(masses[0],2)*(pow(masses[1],4) - 6*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4))) - 2*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],4)*log(pow(masses[0],2)/pow(masses[2],2))))/(2*pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 97: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},58,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},54,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},47,mubarsq); }
            else { return (2*(-3*pow(masses[1],4) + pow(masses[0],2)*pow(masses[2],2) + 2*pow(masses[1],2)*pow(masses[2],2))*log(pow(masses[0],2)/pow(masses[1],2)) + ((pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(2*pow(masses[0],6) - pow(masses[1],4)*pow(masses[2],2) + pow(masses[1],2)*pow(masses[2],4) + pow(masses[0],4)*(5*pow(masses[1],2) - 9*pow(masses[2],2)) - pow(masses[0],2)*(pow(masses[1],4) + 2*pow(masses[1],2)*pow(masses[2],2) - 5*pow(masses[2],4))) - 2*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[2],2),3)))/(2*pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[1],2) - pow(masses[2],2),2)); }
            break;
        case 98: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},55,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},55,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},51,mubarsq); }
            else { return (6*pow(masses[1],8)*pow(pow(masses[0],2) - pow(masses[2],2),4)*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(2*pow(masses[0],8) + 26*pow(masses[1],4)*pow(masses[2],4) - 7*pow(masses[0],6)*(pow(masses[1],2) + pow(masses[2],2)) - 31*pow(masses[0],2)*pow(masses[1],2)*pow(masses[2],2)*(pow(masses[1],2) + pow(masses[2],2)) + pow(masses[0],4)*(11*pow(masses[1],4) + 26*pow(masses[1],2)*pow(masses[2],2) + 11*pow(masses[2],4))) + 6*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],8)*log(pow(masses[0],2)/pow(masses[2],2))))/(6*pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[0],2) - pow(masses[2],2),4)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 99: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},56,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},56,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},52,mubarsq); }
            else { return (6*pow(masses[1],6)*pow(pow(masses[0],2) - pow(masses[2],2),4)*log(pow(masses[0],2)/pow(masses[1],2)) + (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[0],8) - 11*pow(masses[1],4)*pow(masses[2],4) - 5*pow(masses[0],6)*(pow(masses[1],2) + pow(masses[2],2)) + 7*pow(masses[0],2)*pow(masses[1],2)*pow(masses[2],2)*(pow(masses[1],2) + pow(masses[2],2)) - 2*pow(masses[0],4)*(pow(masses[1],4) - 5*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4))) - 6*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],6)*log(pow(masses[0],2)/pow(masses[2],2))))/(6*pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[0],2) - pow(masses[2],2),4)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 100: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},57,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},57,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},53,mubarsq); }
            else { return (6*pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)) - ((pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(2*pow(masses[0],8) + 2*pow(masses[1],4)*pow(masses[2],4) + 5*pow(masses[0],6)*(pow(masses[1],2) + pow(masses[2],2)) + 5*pow(masses[0],2)*pow(masses[1],2)*pow(masses[2],2)*(pow(masses[1],2) + pow(masses[2],2)) - pow(masses[0],4)*(pow(masses[1],4) + 22*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4))) + 6*pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],4)*log(pow(masses[0],2)/pow(masses[2],2))))/(pow(masses[0],2)*pow(pow(masses[0],2) - pow(masses[2],2),4)))/(6*pow(pow(masses[0],2) - pow(masses[1],2),4)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 101: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2]},58,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},58,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1]},54,mubarsq); }
            else { return (6*pow(masses[1],2)*pow(pow(masses[0],3) - masses[0]*pow(masses[2],2),4)*log(pow(masses[0],2)/pow(masses[1],2)) - (pow(masses[0],2) - pow(masses[1],2))*((pow(masses[0],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[2],2))*(11*pow(masses[0],8) - pow(masses[1],4)*pow(masses[2],4) - 7*pow(masses[0],6)*(pow(masses[1],2) + pow(masses[2],2)) + 5*pow(masses[0],2)*pow(masses[1],2)*pow(masses[2],2)*(pow(masses[1],2) + pow(masses[2],2)) + 2*pow(masses[0],4)*(pow(masses[1],4) - 5*pow(masses[1],2)*pow(masses[2],2) + pow(masses[2],4))) + 6*pow(masses[0],4)*pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(masses[2],2)*log(pow(masses[0],2)/pow(masses[2],2))))/(6*pow(masses[0],4)*pow(pow(masses[0],2) - pow(masses[1],2),4)*pow(pow(masses[0],2) - pow(masses[2],2),4)*(pow(masses[1],2) - pow(masses[2],2))); }
            break;
        case 102: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},76,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},76,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},76,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},66,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},66,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},61,mubarsq); }
            else { return ((pow(masses[1],4)*log(pow(masses[0],2)/pow(masses[1],2)))/(-pow(masses[0],2) + pow(masses[1],2)) + ((pow(masses[2],4)*(pow(masses[1],2) - pow(masses[3],2))*log(pow(masses[0],2)/pow(masses[2],2)))/(pow(masses[0],2) - pow(masses[2],2)) + ((pow(masses[1],2) - pow(masses[2],2))*pow(masses[3],4)*log(pow(masses[0],2)/pow(masses[3],2)))/(-pow(masses[0],2) + pow(masses[3],2)))/(pow(masses[2],2) - pow(masses[3],2)))/((pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))); }
            break;
        case 103: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},77,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},77,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},77,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},67,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},67,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},62,mubarsq); }
            else { return ((pow(masses[1],2)*log(pow(masses[0],2)/pow(masses[1],2)))/(-pow(masses[0],2) + pow(masses[1],2)) + ((pow(masses[2],2)*(pow(masses[1],2) - pow(masses[3],2))*log(pow(masses[0],2)/pow(masses[2],2)))/(pow(masses[0],2) - pow(masses[2],2)) + ((pow(masses[1],2) - pow(masses[2],2))*pow(masses[3],2)*log(pow(masses[0],2)/pow(masses[3],2)))/(-pow(masses[0],2) + pow(masses[3],2)))/(pow(masses[2],2) - pow(masses[3],2)))/((pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))); }
            break;
        case 104: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},79,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},79,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},90,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},68,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},71,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},63,mubarsq); }
            else { return pow(masses[0],4)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[0],2) - pow(masses[3],2),2)) + (pow(masses[1],4)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[3],2),2)) + (pow(masses[2],4)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*pow(pow(masses[2],2) - pow(masses[3],2),2)) - (pow(masses[3],4)*(-2*pow(masses[2],2)*pow(masses[3],2) + 3*pow(masses[3],4) + pow(masses[1],2)*(pow(masses[2],2) - 2*pow(masses[3],2)) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2) - 2*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)) + (pow(masses[3],2)*(1 + 2*log(pow(masses[0],2)/pow(masses[3],2))))/((-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 105: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},80,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},80,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},91,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},69,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},72,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},64,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[0],2) - pow(masses[3],2),2)) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[3],2),2)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*pow(pow(masses[2],2) - pow(masses[3],2),2)) + log(pow(masses[0],2)/pow(masses[3],2))/((-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))) - (pow(masses[3],2)*(-2*pow(masses[2],2)*pow(masses[3],2) + 3*pow(masses[3],4) + pow(masses[1],2)*(pow(masses[2],2) - 2*pow(masses[3],2)) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2) - 2*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)); }
            break;
        case 106: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},81,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},81,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},101,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},70,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},74,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},65,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[0],2) - pow(masses[3],2),3)) - 1/(2*pow(masses[3],2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[3],2),3)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*pow(pow(masses[2],2) - pow(masses[3],2),3)) - ((-2*pow(masses[2],2)*pow(masses[3],2) + 3*pow(masses[3],4) + pow(masses[1],2)*(pow(masses[2],2) - 2*pow(masses[3],2)) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2) - 2*pow(masses[3],2)))*log(pow(masses[0],2)/pow(masses[3],2)))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)) - (pow(masses[3],2)*(3*pow(masses[2],4)*pow(masses[3],4) - 8*pow(masses[2],2)*pow(masses[3],6) + 6*pow(masses[3],8) + pow(masses[1],4)*(pow(masses[2],4) - 3*pow(masses[2],2)*pow(masses[3],2) + 3*pow(masses[3],4)) + pow(masses[1],2)*(-3*pow(masses[2],4)*pow(masses[3],2) + 9*pow(masses[2],2)*pow(masses[3],4) - 8*pow(masses[3],6)) + pow(masses[0],4)*(pow(masses[1],4) + pow(masses[2],4) - 3*pow(masses[2],2)*pow(masses[3],2) + 3*pow(masses[3],4) + pow(masses[1],2)*(pow(masses[2],2) - 3*pow(masses[3],2))) + pow(masses[0],2)*(-3*pow(masses[2],4)*pow(masses[3],2) + 9*pow(masses[2],2)*pow(masses[3],4) - 8*pow(masses[3],6) + pow(masses[1],4)*(pow(masses[2],2) - 3*pow(masses[3],2)) + pow(masses[1],2)*pow(pow(masses[2],2) - 3*pow(masses[3],2),2)))*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),3)*pow(-pow(masses[1],2) + pow(masses[3],2),3)*pow(-pow(masses[2],2) + pow(masses[3],2),3)); }
            break;
        case 107: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},83,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},90,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},79,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},71,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},68,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},63,mubarsq); }
            else { return pow(masses[0],4)/((pow(masses[0],2) - pow(masses[1],2))*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[0],2) - pow(masses[3],2))) + (pow(masses[1],4)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*pow(pow(masses[1],2) - pow(masses[2],2),2)*(pow(masses[1],2) - pow(masses[3],2))) - (pow(masses[2],4)*(3*pow(masses[2],4) - 2*pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(-2*pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(pow(masses[1],2) - 2*pow(masses[2],2) + pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)) + (pow(masses[2],2)*(1 + 2*log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],4)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[2],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))); }
            break;
        case 108: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},84,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},91,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},80,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},72,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},69,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},64,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[0],2) - pow(masses[3],2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*pow(pow(masses[1],2) - pow(masses[2],2),2)*(pow(masses[1],2) - pow(masses[3],2))) + log(pow(masses[0],2)/pow(masses[2],2))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) - (pow(masses[2],2)*(3*pow(masses[2],4) - 2*pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(-2*pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(pow(masses[1],2) - 2*pow(masses[2],2) + pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[2],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))); }
            break;
        case 109: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},86,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},93,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},93,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},73,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},73,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},65,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),2)) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)) + log(pow(masses[0],2)/pow(masses[2],2))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*pow(pow(masses[2],2) - pow(masses[3],2),2)) - (pow(masses[2],2)*(4*pow(masses[2],4) - 2*pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(-3*pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(2*pow(masses[1],2) - 3*pow(masses[2],2) + pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),3)) + log(pow(masses[0],2)/pow(masses[3],2))/(pow(pow(masses[2],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))) - (pow(masses[3],2)*(-2*pow(masses[2],2)*pow(masses[3],2) + 4*pow(masses[3],4) + pow(masses[1],2)*(pow(masses[2],2) - 3*pow(masses[3],2)) + pow(masses[0],2)*(2*pow(masses[1],2) + pow(masses[2],2) - 3*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)*pow(-pow(masses[2],2) + pow(masses[3],2),3)); }
            break;
        case 110: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},88,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},101,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},81,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},74,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},70,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},65,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*pow(pow(masses[0],2) - pow(masses[2],2),3)*(pow(masses[0],2) - pow(masses[3],2))) - 1/(2*pow(masses[2],2)*(-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*pow(pow(masses[1],2) - pow(masses[2],2),3)*(pow(masses[1],2) - pow(masses[3],2))) - ((3*pow(masses[2],4) - 2*pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(-2*pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(pow(masses[1],2) - 2*pow(masses[2],2) + pow(masses[3],2)))*log(pow(masses[0],2)/pow(masses[2],2)))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)) - (pow(masses[2],2)*(6*pow(masses[2],8) - 8*pow(masses[2],6)*pow(masses[3],2) + 3*pow(masses[2],4)*pow(masses[3],4) + pow(masses[1],4)*(3*pow(masses[2],4) - 3*pow(masses[2],2)*pow(masses[3],2) + pow(masses[3],4)) + pow(masses[1],2)*(-8*pow(masses[2],6) + 9*pow(masses[2],4)*pow(masses[3],2) - 3*pow(masses[2],2)*pow(masses[3],4)) + pow(masses[0],4)*(pow(masses[1],4) + 3*pow(masses[2],4) - 3*pow(masses[2],2)*pow(masses[3],2) + pow(masses[3],4) + pow(masses[1],2)*(-3*pow(masses[2],2) + pow(masses[3],2))) + pow(masses[0],2)*(-8*pow(masses[2],6) + 9*pow(masses[2],4)*pow(masses[3],2) - 3*pow(masses[2],2)*pow(masses[3],4) + pow(masses[1],4)*(-3*pow(masses[2],2) + pow(masses[3],2)) + pow(masses[1],2)*pow(-3*pow(masses[2],2) + pow(masses[3],2),2)))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(-pow(masses[1],2) + pow(masses[2],2),3)*pow(pow(masses[2],2) - pow(masses[3],2),3)) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/((-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*pow(-pow(masses[2],2) + pow(masses[3],2),3)); }
            break;
        case 111: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},90,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},83,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},83,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},71,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},71,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},68,mubarsq); }
            else { return pow(masses[0],4)/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))) - (pow(masses[1],4)*(3*pow(masses[1],4) + pow(masses[2],2)*pow(masses[3],2) - 2*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(-2*pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)) + (pow(masses[1],2)*(1 + 2*log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + (pow(masses[2],4)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[1],2) - pow(masses[2],2),2)*(-pow(masses[0],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],4)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[1],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 112: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},91,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},84,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},84,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},72,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},72,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},69,mubarsq); }
            else { return pow(masses[0],2)/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))) + log(pow(masses[0],2)/pow(masses[1],2))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) - (pow(masses[1],2)*(3*pow(masses[1],4) + pow(masses[2],2)*pow(masses[3],2) - 2*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(-2*pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[1],2) - pow(masses[2],2),2)*(-pow(masses[0],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[1],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 113: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},93,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},86,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},97,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},73,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},74,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},70,mubarsq); }
            else { return pow(masses[0],2)/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[0],2) - pow(masses[3],2),2)) + log(pow(masses[0],2)/pow(masses[1],2))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[3],2),2)) - (pow(masses[1],2)*(4*pow(masses[1],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[0],2)*(-3*pow(masses[1],2) + 2*pow(masses[2],2) + pow(masses[3],2)) - pow(masses[1],2)*(3*pow(masses[2],2) + 2*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),3)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[1],2) - pow(masses[2],2),2)*(-pow(masses[0],2) + pow(masses[2],2))*pow(pow(masses[2],2) - pow(masses[3],2),2)) + log(pow(masses[0],2)/pow(masses[3],2))/(pow(pow(masses[1],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))) - (pow(masses[3],2)*(-3*pow(masses[2],2)*pow(masses[3],2) + 4*pow(masses[3],4) + pow(masses[0],2)*(pow(masses[1],2) + 2*pow(masses[2],2) - 3*pow(masses[3],2)) + pow(masses[1],2)*(pow(masses[2],2) - 2*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)*pow(-pow(masses[1],2) + pow(masses[3],2),3)); }
            break;
        case 114: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},97,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},97,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},86,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},74,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},73,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},70,mubarsq); }
            else { return pow(masses[0],2)/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[0],2) - pow(masses[3],2))) + log(pow(masses[0],2)/pow(masses[1],2))/((-pow(masses[0],2) + pow(masses[1],2))*pow(pow(masses[1],2) - pow(masses[2],2),2)*(pow(masses[1],2) - pow(masses[3],2))) - (pow(masses[1],2)*(4*pow(masses[1],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[0],2)*(-3*pow(masses[1],2) + pow(masses[2],2) + 2*pow(masses[3],2)) - pow(masses[1],2)*(2*pow(masses[2],2) + 3*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[3],2),2)) + log(pow(masses[0],2)/pow(masses[2],2))/(pow(pow(masses[1],2) - pow(masses[2],2),2)*(-pow(masses[0],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) - (pow(masses[2],2)*(4*pow(masses[2],4) - 3*pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(-2*pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(pow(masses[1],2) - 3*pow(masses[2],2) + 2*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(-pow(masses[1],2) + pow(masses[2],2),3)*pow(pow(masses[2],2) - pow(masses[3],2),2)) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[1],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))); }
            break;
        case 115: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},101,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},88,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},88,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},74,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},74,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},73,mubarsq); }
            else { return pow(masses[0],2)/(pow(pow(masses[0],2) - pow(masses[1],2),3)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))) - 1/(2*pow(masses[1],2)*(-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + ((-3*pow(masses[1],4) - pow(masses[2],2)*pow(masses[3],2) + pow(masses[0],2)*(2*pow(masses[1],2) - pow(masses[2],2) - pow(masses[3],2)) + 2*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)))*log(pow(masses[0],2)/pow(masses[1],2)))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)) - (pow(masses[1],2)*(6*pow(masses[1],8) + pow(masses[2],4)*pow(masses[3],4) - 8*pow(masses[1],6)*(pow(masses[2],2) + pow(masses[3],2)) - 3*pow(masses[1],2)*pow(masses[2],2)*pow(masses[3],2)*(pow(masses[2],2) + pow(masses[3],2)) + 3*pow(masses[1],4)*(pow(masses[2],4) + 3*pow(masses[2],2)*pow(masses[3],2) + pow(masses[3],4)) + pow(masses[0],4)*(3*pow(masses[1],4) + pow(masses[2],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[3],4) - 3*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2))) + pow(masses[0],2)*(-8*pow(masses[1],6) + 9*pow(masses[1],4)*(pow(masses[2],2) + pow(masses[3],2)) + pow(masses[2],2)*pow(masses[3],2)*(pow(masses[2],2) + pow(masses[3],2)) - 3*pow(masses[1],2)*pow(pow(masses[2],2) + pow(masses[3],2),2)))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[3],2),3)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*pow(-pow(masses[1],2) + pow(masses[2],2),3)*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/((-pow(masses[0],2) + pow(masses[3],2))*pow(-pow(masses[1],2) + pow(masses[3],2),3)*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 116: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},89,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},89,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},89,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},82,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},82,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},78,mubarsq); }
            else { return (2*pow(masses[0],4))/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))) - (pow(masses[0],6)*(3*pow(masses[0],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) - 2*pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),2)) + (pow(masses[1],6)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + (pow(masses[2],6)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],6)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 117: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},90,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},90,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},90,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},83,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},83,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},79,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))) - (pow(masses[0],4)*(3*pow(masses[0],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) - 2*pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),2)) + (pow(masses[1],4)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + (pow(masses[2],4)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],4)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 118: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},91,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},91,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},91,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},84,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},84,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},80,mubarsq); }
            else { return -((pow(masses[0],2)*(3*pow(masses[0],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) - 2*pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 119: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},93,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},93,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},101,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},86,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},88,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},81,mubarsq); }
            else { return -((pow(masses[0],2)*(4*pow(masses[0],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(2*pow(masses[2],2) + pow(masses[3],2)) - pow(masses[0],2)*(3*pow(masses[1],2) + 3*pow(masses[2],2) + 2*pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),3))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[1],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[3],2),2)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*(-pow(masses[1],2) + pow(masses[2],2))*pow(pow(masses[2],2) - pow(masses[3],2),2)) + log(pow(masses[0],2)/pow(masses[3],2))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))) + (pow(masses[3],2)*(-3*pow(masses[2],2)*pow(masses[3],2) + 4*pow(masses[3],4) + pow(masses[1],2)*(2*pow(masses[2],2) - 3*pow(masses[3],2)) + pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2) - 2*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),3)*pow(pow(masses[1],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)); }
            break;
        case 120: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},97,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},101,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},93,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},88,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},86,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},81,mubarsq); }
            else { return -((pow(masses[0],2)*(4*pow(masses[0],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(pow(masses[2],2) + 2*pow(masses[3],2)) - pow(masses[0],2)*(3*pow(masses[1],2) + 2*pow(masses[2],2) + 3*pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[0],2) - pow(masses[3],2),2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*(pow(masses[1],2) - pow(masses[3],2))) + log(pow(masses[0],2)/pow(masses[2],2))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[2],2)*(4*pow(masses[2],4) - 3*pow(masses[2],2)*pow(masses[3],2) + pow(masses[0],2)*(pow(masses[1],2) - 2*pow(masses[2],2) + pow(masses[3],2)) + pow(masses[1],2)*(-3*pow(masses[2],2) + 2*pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)*(-pow(masses[1],2) + pow(masses[3],2))); }
            break;
        case 121: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},100,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},96,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},96,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},87,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},87,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},85,mubarsq); }
            else { return pow(masses[0],2)/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))) - (pow(masses[0],4)*(4*pow(masses[0],4) + 2*pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) - pow(masses[0],2)*(2*pow(masses[1],2) + 3*(pow(masses[2],2) + pow(masses[3],2)))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),2)) + (pow(masses[1],4)*(4*pow(masses[1],4) + 2*pow(masses[2],2)*pow(masses[3],2) - 3*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(-2*pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)) + (pow(masses[1],2)*(1 + 2*log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + (pow(masses[2],4)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],4)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 122: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},101,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},97,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},97,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},88,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},88,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},86,mubarsq); }
            else { return -((pow(masses[0],2)*(4*pow(masses[0],4) + 2*pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) - pow(masses[0],2)*(2*pow(masses[1],2) + 3*(pow(masses[2],2) + pow(masses[3],2)))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),2))) + log(pow(masses[0],2)/pow(masses[1],2))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + (pow(masses[1],2)*(4*pow(masses[1],4) + 2*pow(masses[2],2)*pow(masses[3],2) - 3*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) + pow(masses[0],2)*(-2*pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2)))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 123: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},100,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},100,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},100,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},96,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},96,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},92,mubarsq); }
            else { return -0.5*1/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))) - (pow(masses[0],2)*(3*pow(masses[0],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)) - 2*pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),2)) + (pow(masses[0],4)*(6*pow(masses[0],8) + pow(masses[2],4)*pow(masses[3],4) + pow(masses[1],2)*pow(masses[2],2)*pow(masses[3],2)*(pow(masses[2],2) + pow(masses[3],2)) - 3*pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))*(pow(masses[1],2) + pow(masses[3],2))*(pow(masses[2],2) + pow(masses[3],2)) - 8*pow(masses[0],6)*(pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2)) + pow(masses[1],4)*(pow(masses[2],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[3],4)) + 3*pow(masses[0],4)*(pow(masses[1],4) + pow(masses[2],4) + 3*pow(masses[2],2)*pow(masses[3],2) + pow(masses[3],4) + 3*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[0],2) - pow(masses[3],2),3)) + (pow(masses[1],4)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(-pow(masses[0],2) + pow(masses[1],2),3)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + (pow(masses[2],4)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(-pow(masses[0],2) + pow(masses[2],2),3)*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],4)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(-pow(masses[0],2) + pow(masses[3],2),3)*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 124: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3]},101,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},101,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},101,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3]},97,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},97,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2]},93,mubarsq); }
            else { return -0.5*1/(pow(masses[0],2)*(pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))) + (pow(masses[0],2)*(6*pow(masses[0],8) + pow(masses[2],4)*pow(masses[3],4) + pow(masses[1],2)*pow(masses[2],2)*pow(masses[3],2)*(pow(masses[2],2) + pow(masses[3],2)) - 3*pow(masses[0],2)*(pow(masses[1],2) + pow(masses[2],2))*(pow(masses[1],2) + pow(masses[3],2))*(pow(masses[2],2) + pow(masses[3],2)) - 8*pow(masses[0],6)*(pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2)) + pow(masses[1],4)*(pow(masses[2],4) + pow(masses[2],2)*pow(masses[3],2) + pow(masses[3],4)) + 3*pow(masses[0],4)*(pow(masses[1],4) + pow(masses[2],4) + 3*pow(masses[2],2)*pow(masses[3],2) + pow(masses[3],4) + 3*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2)))))/(pow(pow(masses[0],2) - pow(masses[1],2),3)*pow(pow(masses[0],2) - pow(masses[2],2),3)*pow(pow(masses[0],2) - pow(masses[3],2),3)) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(-pow(masses[0],2) + pow(masses[1],2),3)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(-pow(masses[0],2) + pow(masses[2],2),3)*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(-pow(masses[0],2) + pow(masses[3],2),3)*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))); }
            break;
        case 125: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3],masses[4]},117,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},117,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},117,mubarsq); }
            else if (abs(masses[0] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},117,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},111,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},111,mubarsq); }
            else if (abs(masses[1] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},111,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},107,mubarsq); }
            else if (abs(masses[2] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},107,mubarsq); }
            else if (abs(masses[3] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},104,mubarsq); }
            else { return pow(masses[0],4)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))*(pow(masses[0],2) - pow(masses[4],2))) + (pow(masses[1],4)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))*(pow(masses[1],2) - pow(masses[4],2))) + (pow(masses[2],4)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))*(pow(masses[2],2) - pow(masses[4],2))) + (pow(masses[3],4)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/((-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))*(pow(masses[3],2) - pow(masses[4],2))) + (pow(masses[4],4)*(1 + log(pow(masses[0],2)/pow(masses[4],2))))/((-pow(masses[0],2) + pow(masses[4],2))*(-pow(masses[1],2) + pow(masses[4],2))*(-pow(masses[2],2) + pow(masses[4],2))*(-pow(masses[3],2) + pow(masses[4],2))); }
            break;
        case 126: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3],masses[4]},118,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},118,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},118,mubarsq); }
            else if (abs(masses[0] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},118,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},112,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},112,mubarsq); }
            else if (abs(masses[1] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},112,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},108,mubarsq); }
            else if (abs(masses[2] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},108,mubarsq); }
            else if (abs(masses[3] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},105,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))*(pow(masses[0],2) - pow(masses[4],2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))*(pow(masses[1],2) - pow(masses[4],2))) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))*(pow(masses[2],2) - pow(masses[4],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/((-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))*(pow(masses[3],2) - pow(masses[4],2))) + (pow(masses[4],2)*(1 + log(pow(masses[0],2)/pow(masses[4],2))))/((-pow(masses[0],2) + pow(masses[4],2))*(-pow(masses[1],2) + pow(masses[4],2))*(-pow(masses[2],2) + pow(masses[4],2))*(-pow(masses[3],2) + pow(masses[4],2))); }
            break;
        case 127: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3],masses[4]},119,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},119,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},119,mubarsq); }
            else if (abs(masses[0] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},124,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},113,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},113,mubarsq); }
            else if (abs(masses[1] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},115,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},109,mubarsq); }
            else if (abs(masses[2] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},110,mubarsq); }
            else if (abs(masses[3] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},106,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))*pow(pow(masses[0],2) - pow(masses[4],2),2)) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))*pow(pow(masses[1],2) - pow(masses[4],2),2)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))*pow(pow(masses[2],2) - pow(masses[4],2),2)) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/((-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))*pow(pow(masses[3],2) - pow(masses[4],2),2)) + log(pow(masses[0],2)/pow(masses[4],2))/((-pow(masses[0],2) + pow(masses[4],2))*(-pow(masses[1],2) + pow(masses[4],2))*(-pow(masses[2],2) + pow(masses[4],2))*(-pow(masses[3],2) + pow(masses[4],2))) - (pow(masses[4],2)*(2*pow(masses[2],2)*pow(masses[3],2)*pow(masses[4],2) - 3*pow(masses[2],2)*pow(masses[4],4) - 3*pow(masses[3],2)*pow(masses[4],4) + 4*pow(masses[4],6) + pow(masses[1],2)*(2*pow(masses[3],2)*pow(masses[4],2) - 3*pow(masses[4],4) - pow(masses[2],2)*(pow(masses[3],2) - 2*pow(masses[4],2))) - pow(masses[0],2)*(-2*pow(masses[3],2)*pow(masses[4],2) + 3*pow(masses[4],4) + pow(masses[2],2)*(pow(masses[3],2) - 2*pow(masses[4],2)) + pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2) - 2*pow(masses[4],2))))*(1 + log(pow(masses[0],2)/pow(masses[4],2))))/(pow(pow(masses[0],2) - pow(masses[4],2),2)*pow(pow(masses[1],2) - pow(masses[4],2),2)*pow(pow(masses[2],2) - pow(masses[4],2),2)*pow(pow(masses[3],2) - pow(masses[4],2),2)); }
            break;
        case 128: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3],masses[4]},120,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},120,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},124,mubarsq); }
            else if (abs(masses[0] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},119,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},114,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},115,mubarsq); }
            else if (abs(masses[1] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},113,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},110,mubarsq); }
            else if (abs(masses[2] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},109,mubarsq); }
            else if (abs(masses[3] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},106,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*pow(pow(masses[0],2) - pow(masses[3],2),2)*(pow(masses[0],2) - pow(masses[4],2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*pow(pow(masses[1],2) - pow(masses[3],2),2)*(pow(masses[1],2) - pow(masses[4],2))) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*pow(pow(masses[2],2) - pow(masses[3],2),2)*(pow(masses[2],2) - pow(masses[4],2))) + log(pow(masses[0],2)/pow(masses[3],2))/((-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))*(pow(masses[3],2) - pow(masses[4],2))) - (pow(masses[3],2)*(-3*pow(masses[2],2)*pow(masses[3],4) + 4*pow(masses[3],6) + 2*pow(masses[2],2)*pow(masses[3],2)*pow(masses[4],2) - 3*pow(masses[3],4)*pow(masses[4],2) + pow(masses[1],2)*(-3*pow(masses[3],4) + 2*pow(masses[3],2)*pow(masses[4],2) + pow(masses[2],2)*(2*pow(masses[3],2) - pow(masses[4],2))) - pow(masses[0],2)*(3*pow(masses[3],4) - 2*pow(masses[3],2)*pow(masses[4],2) + pow(masses[2],2)*(-2*pow(masses[3],2) + pow(masses[4],2)) + pow(masses[1],2)*(pow(masses[2],2) - 2*pow(masses[3],2) + pow(masses[4],2))))*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)*pow(pow(masses[3],2) - pow(masses[4],2),2)) + (pow(masses[4],2)*(1 + log(pow(masses[0],2)/pow(masses[4],2))))/(pow(pow(masses[3],2) - pow(masses[4],2),2)*(-pow(masses[0],2) + pow(masses[4],2))*(-pow(masses[1],2) + pow(masses[4],2))*(-pow(masses[2],2) + pow(masses[4],2))); }
            break;
        case 129: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3],masses[4]},122,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},124,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},120,mubarsq); }
            else if (abs(masses[0] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},120,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},115,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},114,mubarsq); }
            else if (abs(masses[1] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},114,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},110,mubarsq); }
            else if (abs(masses[2] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},110,mubarsq); }
            else if (abs(masses[3] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},109,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*pow(pow(masses[0],2) - pow(masses[2],2),2)*(pow(masses[0],2) - pow(masses[3],2))*(pow(masses[0],2) - pow(masses[4],2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*pow(pow(masses[1],2) - pow(masses[2],2),2)*(pow(masses[1],2) - pow(masses[3],2))*(pow(masses[1],2) - pow(masses[4],2))) + log(pow(masses[0],2)/pow(masses[2],2))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))*(pow(masses[2],2) - pow(masses[4],2))) - (pow(masses[2],2)*(4*pow(masses[2],6) - 3*pow(masses[2],4)*pow(masses[3],2) - 3*pow(masses[2],4)*pow(masses[4],2) + 2*pow(masses[2],2)*pow(masses[3],2)*pow(masses[4],2) + pow(masses[1],2)*(-3*pow(masses[2],4) - pow(masses[3],2)*pow(masses[4],2) + 2*pow(masses[2],2)*(pow(masses[3],2) + pow(masses[4],2))) + pow(masses[0],2)*(-3*pow(masses[2],4) - pow(masses[3],2)*pow(masses[4],2) + pow(masses[1],2)*(2*pow(masses[2],2) - pow(masses[3],2) - pow(masses[4],2)) + 2*pow(masses[2],2)*(pow(masses[3],2) + pow(masses[4],2))))*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[2],2) - pow(masses[3],2),2)*pow(pow(masses[2],2) - pow(masses[4],2),2)) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[2],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(pow(masses[3],2) - pow(masses[4],2))) + (pow(masses[4],2)*(1 + log(pow(masses[0],2)/pow(masses[4],2))))/(pow(pow(masses[2],2) - pow(masses[4],2),2)*(-pow(masses[0],2) + pow(masses[4],2))*(-pow(masses[1],2) + pow(masses[4],2))*(-pow(masses[3],2) + pow(masses[4],2))); }
            break;
        case 130: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3],masses[4]},124,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},122,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},122,mubarsq); }
            else if (abs(masses[0] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},122,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},115,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},115,mubarsq); }
            else if (abs(masses[1] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},115,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},114,mubarsq); }
            else if (abs(masses[2] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},114,mubarsq); }
            else if (abs(masses[3] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},113,mubarsq); }
            else { return pow(masses[0],2)/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))*(pow(masses[0],2) - pow(masses[4],2))) + log(pow(masses[0],2)/pow(masses[1],2))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))*(pow(masses[1],2) - pow(masses[4],2))) - (pow(masses[1],2)*(4*pow(masses[1],6) - pow(masses[2],2)*pow(masses[3],2)*pow(masses[4],2) - 3*pow(masses[1],4)*(pow(masses[2],2) + pow(masses[3],2) + pow(masses[4],2)) + 2*pow(masses[1],2)*(pow(masses[3],2)*pow(masses[4],2) + pow(masses[2],2)*(pow(masses[3],2) + pow(masses[4],2))) - pow(masses[0],2)*(3*pow(masses[1],4) + pow(masses[3],2)*pow(masses[4],2) + pow(masses[2],2)*(pow(masses[3],2) + pow(masses[4],2)) - 2*pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2) + pow(masses[4],2))))*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[1],2) - pow(masses[2],2),2)*pow(pow(masses[1],2) - pow(masses[3],2),2)*pow(pow(masses[1],2) - pow(masses[4],2),2)) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[1],2) - pow(masses[2],2),2)*(-pow(masses[0],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))*(pow(masses[2],2) - pow(masses[4],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[1],2) - pow(masses[3],2),2)*(-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))*(pow(masses[3],2) - pow(masses[4],2))) + (pow(masses[4],2)*(1 + log(pow(masses[0],2)/pow(masses[4],2))))/(pow(pow(masses[1],2) - pow(masses[4],2),2)*(-pow(masses[0],2) + pow(masses[4],2))*(-pow(masses[2],2) + pow(masses[4],2))*(-pow(masses[3],2) + pow(masses[4],2))); }
            break;
        case 131: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3],masses[4]},124,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},124,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},124,mubarsq); }
            else if (abs(masses[0] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},124,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4]},122,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},122,mubarsq); }
            else if (abs(masses[1] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},122,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4]},120,mubarsq); }
            else if (abs(masses[2] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},120,mubarsq); }
            else if (abs(masses[3] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3]},119,mubarsq); }
            else { return (pow(masses[0],2)*(-4*pow(masses[0],6) + pow(masses[2],2)*pow(masses[3],2)*pow(masses[4],2) + 3*pow(masses[0],4)*(pow(masses[1],2) + pow(masses[2],2) + pow(masses[3],2) + pow(masses[4],2)) + pow(masses[1],2)*(pow(masses[3],2)*pow(masses[4],2) + pow(masses[2],2)*(pow(masses[3],2) + pow(masses[4],2))) - 2*pow(masses[0],2)*(pow(masses[3],2)*pow(masses[4],2) + pow(masses[2],2)*(pow(masses[3],2) + pow(masses[4],2)) + pow(masses[1],2)*(pow(masses[2],2) + pow(masses[3],2) + pow(masses[4],2)))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*pow(pow(masses[0],2) - pow(masses[2],2),2)*pow(pow(masses[0],2) - pow(masses[3],2),2)*pow(pow(masses[0],2) - pow(masses[4],2),2)) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/(pow(pow(masses[0],2) - pow(masses[1],2),2)*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))*(pow(masses[1],2) - pow(masses[4],2))) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/(pow(pow(masses[0],2) - pow(masses[2],2),2)*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))*(pow(masses[2],2) - pow(masses[4],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/(pow(pow(masses[0],2) - pow(masses[3],2),2)*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))*(pow(masses[3],2) - pow(masses[4],2))) + (pow(masses[4],2)*(1 + log(pow(masses[0],2)/pow(masses[4],2))))/(pow(pow(masses[0],2) - pow(masses[4],2),2)*(-pow(masses[1],2) + pow(masses[4],2))*(-pow(masses[2],2) + pow(masses[4],2))*(-pow(masses[3],2) + pow(masses[4],2))); }
            break;
        case 132: 
            if (abs(masses[0] - masses[1]) < 1e-6) { return LF({masses[0],masses[2],masses[3],masses[4],masses[5]},131,mubarsq); }
            else if (abs(masses[0] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4],masses[5]},131,mubarsq); }
            else if (abs(masses[0] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4],masses[5]},131,mubarsq); }
            else if (abs(masses[0] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[5]},131,mubarsq); }
            else if (abs(masses[0] - masses[5]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[4]},131,mubarsq); }
            else if (abs(masses[1] - masses[2]) < 1e-6) { return LF({masses[0],masses[1],masses[3],masses[4],masses[5]},130,mubarsq); }
            else if (abs(masses[1] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4],masses[5]},130,mubarsq); }
            else if (abs(masses[1] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[5]},130,mubarsq); }
            else if (abs(masses[1] - masses[5]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[4]},130,mubarsq); }
            else if (abs(masses[2] - masses[3]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[4],masses[5]},129,mubarsq); }
            else if (abs(masses[2] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[5]},129,mubarsq); }
            else if (abs(masses[2] - masses[5]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[4]},129,mubarsq); }
            else if (abs(masses[3] - masses[4]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[5]},128,mubarsq); }
            else if (abs(masses[3] - masses[5]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[4]},128,mubarsq); }
            else if (abs(masses[4] - masses[5]) < 1e-6) { return LF({masses[0],masses[1],masses[2],masses[3],masses[4]},127,mubarsq); }
            else { return pow(masses[0],2)/((pow(masses[0],2) - pow(masses[1],2))*(pow(masses[0],2) - pow(masses[2],2))*(pow(masses[0],2) - pow(masses[3],2))*(pow(masses[0],2) - pow(masses[4],2))*(pow(masses[0],2) - pow(masses[5],2))) + (pow(masses[1],2)*(1 + log(pow(masses[0],2)/pow(masses[1],2))))/((-pow(masses[0],2) + pow(masses[1],2))*(pow(masses[1],2) - pow(masses[2],2))*(pow(masses[1],2) - pow(masses[3],2))*(pow(masses[1],2) - pow(masses[4],2))*(pow(masses[1],2) - pow(masses[5],2))) + (pow(masses[2],2)*(1 + log(pow(masses[0],2)/pow(masses[2],2))))/((-pow(masses[0],2) + pow(masses[2],2))*(-pow(masses[1],2) + pow(masses[2],2))*(pow(masses[2],2) - pow(masses[3],2))*(pow(masses[2],2) - pow(masses[4],2))*(pow(masses[2],2) - pow(masses[5],2))) + (pow(masses[3],2)*(1 + log(pow(masses[0],2)/pow(masses[3],2))))/((-pow(masses[0],2) + pow(masses[3],2))*(-pow(masses[1],2) + pow(masses[3],2))*(-pow(masses[2],2) + pow(masses[3],2))*(pow(masses[3],2) - pow(masses[4],2))*(pow(masses[3],2) - pow(masses[5],2))) + (pow(masses[4],2)*(1 + log(pow(masses[0],2)/pow(masses[4],2))))/((-pow(masses[0],2) + pow(masses[4],2))*(-pow(masses[1],2) + pow(masses[4],2))*(-pow(masses[2],2) + pow(masses[4],2))*(-pow(masses[3],2) + pow(masses[4],2))*(pow(masses[4],2) - pow(masses[5],2))) + (pow(masses[5],2)*(1 + log(pow(masses[0],2)/pow(masses[5],2))))/((-pow(masses[0],2) + pow(masses[5],2))*(-pow(masses[1],2) + pow(masses[5],2))*(-pow(masses[2],2) + pow(masses[5],2))*(-pow(masses[3],2) + pow(masses[5],2))*(-pow(masses[4],2) + pow(masses[5],2))); }
            break;
    }
    return 0.0;
}
