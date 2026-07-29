#include "LF_helper.h"
#include "complex_math.h"
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

double rel_diff(std::complex<double> a, std::complex<double> b) {
    return std::abs(a-b) / std::min(std::abs(a),std::abs(b));
}

std::complex<double> LFh1(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)*(1 + log(mubarsq/pow(m[0],2)));
}

std::complex<double> LFh2(const std::vector<std::complex<double> >& m, double mubarsq) {
    return 1 + log(mubarsq/pow(m[0],2));
}

std::complex<double> LFh3(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (1 + log(mubarsq/pow(m[0],2)))/pow(m[0],2);
}

std::complex<double> LFh4(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)*(1 + 2*log(mubarsq/pow(m[0],2)));
}

std::complex<double> LFh5(const std::vector<std::complex<double> >& m, double mubarsq) {
    return log(mubarsq/pow(m[0],2));
}

std::complex<double> LFh6(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -pow(m[0],-2);
}

std::complex<double> LFh7(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5 + log(mubarsq/pow(m[0],2));
}

std::complex<double> LFh8(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*1/pow(m[0],2);
}

std::complex<double> LFh9(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.8333333333333334 + log(mubarsq/pow(m[0],2));
}

std::complex<double> LFh10(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.3333333333333333*1/pow(m[0],2);
}

std::complex<double> LFh11(const std::vector<std::complex<double> >& m, double mubarsq) {
    return 1/(6*pow(m[0],4));
}

std::complex<double> LFh12(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.25*1/pow(m[0],2);
}

std::complex<double> LFh13(const std::vector<std::complex<double> >& m, double mubarsq) {
    return 1/(12*pow(m[0],4));
}

std::complex<double> LFh14(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.08333333333333333*1/pow(m[0],6);
}

std::complex<double> LFh15(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.2*1/pow(m[0],2);
}

std::complex<double> LFh16(const std::vector<std::complex<double> >& m, double mubarsq) {
    return 1/(20*pow(m[0],4));
}

std::complex<double> LFh17(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.03333333333333333*1/pow(m[0],6);
}

std::complex<double> LFh18(const std::vector<std::complex<double> >& m, double mubarsq) {
    return 1/(20*pow(m[0],8));
}

std::complex<double> LFh19(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],4) - pow(m[1],4) + pow(m[0],4)*log(mubarsq/pow(m[0],2)) - pow(m[1],4)*log(mubarsq/pow(m[1],2)))/(pow(m[0],2) - pow(m[1],2));
}

std::complex<double> LFh20(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],2) - pow(m[1],2) + pow(m[0],2)*log(mubarsq/pow(m[0],2)) - pow(m[1],2)*log(mubarsq/pow(m[1],2)))/(pow(m[0],2) - pow(m[1],2));
}

std::complex<double> LFh21(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],4) - pow(m[0],2)*pow(m[1],2) + pow(m[0],4)*log(mubarsq/pow(m[0],2)) + (-2*pow(m[0],2)*pow(m[1],2) + pow(m[1],4))*log(mubarsq/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),2);
}

std::complex<double> LFh22(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],2) - pow(m[1],2) - pow(m[0],2)*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),2);
}

std::complex<double> LFh23(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (3*pow(m[0],4) - 4*pow(m[0],2)*pow(m[1],2) + pow(m[1],4) - 2*pow(m[0],4)*log(pow(m[0],2)/pow(m[1],2)))/(2*pow(pow(m[0],2) - pow(m[1],2),3));
}

std::complex<double> LFh24(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],4) + pow(m[1],4) + 2*pow(m[0],2)*pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/(2*pow(m[1],2)*pow(-pow(m[0],2) + pow(m[1],2),3));
}

std::complex<double> LFh25(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*pow(m[0],6) + 3*pow(m[0],4)*pow(m[1],2) - 6*pow(m[0],2)*pow(m[1],4) + pow(m[1],6) - 6*pow(m[0],4)*pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/(6*pow(m[1],2)*pow(pow(m[0],2) - pow(m[1],2),4));
}

std::complex<double> LFh26(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.16666666666666666*(pow(m[0],6) - 6*pow(m[0],4)*pow(m[1],2) + 3*pow(m[0],2)*pow(m[1],4) + 2*pow(m[1],6) + 6*pow(m[0],2)*pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/pow(-(pow(m[0],2)*m[1]) + pow(m[1],3),4);
}

std::complex<double> LFh27(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],8) + 6*pow(m[0],6)*pow(m[1],2) - 18*pow(m[0],4)*pow(m[1],4) + 10*pow(m[0],2)*pow(m[1],6) + 3*pow(m[1],8) + 12*pow(m[0],2)*pow(m[1],6)*log(pow(m[0],2)/pow(m[1],2)))/(12*pow(m[1],6)*pow(-pow(m[0],2) + pow(m[1],2),5));
}

std::complex<double> LFh28(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-(pow(m[0],2)*pow(m[1],2)) + pow(m[1],4) + (pow(m[0],4) - 2*pow(m[0],2)*pow(m[1],2))*log(mubarsq/pow(m[0],2)) + pow(m[1],4)*log(mubarsq/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),2);
}

std::complex<double> LFh29(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],2) + pow(m[1],2) + pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),2);
}

std::complex<double> LFh30(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[0],4)*pow(m[1],2) + 2*pow(m[0],2)*pow(m[1],4) + (pow(m[0],6) - 3*pow(m[0],4)*pow(m[1],2))*log(mubarsq/pow(m[0],2)) + (3*pow(m[0],2)*pow(m[1],4) - pow(m[1],6))*log(mubarsq/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),3);
}

std::complex<double> LFh31(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],4) + pow(m[1],4) + 2*pow(m[0],2)*pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),3);
}

std::complex<double> LFh32(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[0],2) + 2*pow(m[1],2) + (pow(m[0],2) + pow(m[1],2))*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),3);
}

std::complex<double> LFh33(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*(2*pow(m[0],6) + 3*pow(m[0],4)*pow(m[1],2) - 6*pow(m[0],2)*pow(m[1],4) + pow(m[1],6) - 6*pow(m[0],4)*pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),4);
}

std::complex<double> LFh34(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-5*pow(m[0],4) + 4*pow(m[0],2)*pow(m[1],2) + pow(m[1],4) + 2*(pow(m[0],4) + 2*pow(m[0],2)*pow(m[1],2))*log(pow(m[0],2)/pow(m[1],2)))/(2*pow(pow(m[0],2) - pow(m[1],2),4));
}

std::complex<double> LFh35(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*(pow(m[0],4) + 4*pow(m[0],2)*pow(m[1],2) - 5*pow(m[1],4) - 2*(2*pow(m[0],2)*pow(m[1],2) + pow(m[1],4))*log(pow(m[0],2)/pow(m[1],2)))/(pow(m[1],2)*pow(pow(m[0],2) - pow(m[1],2),4));
}

std::complex<double> LFh36(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],6) + 9*pow(m[0],4)*pow(m[1],2) - 9*pow(m[0],2)*pow(m[1],4) - pow(m[1],6) - 6*pow(m[0],2)*pow(m[1],2)*(pow(m[0],2) + pow(m[1],2))*log(pow(m[0],2)/pow(m[1],2)))/(3*pow(m[1],2)*pow(-pow(m[0],2) + pow(m[1],2),5));
}

std::complex<double> LFh37(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.16666666666666666*(pow(m[0],6) - 9*pow(m[0],4)*pow(m[1],2) - 9*pow(m[0],2)*pow(m[1],4) + 17*pow(m[1],6) + 6*(3*pow(m[0],2)*pow(m[1],4) + pow(m[1],6))*log(pow(m[0],2)/pow(m[1],2)))/(pow(m[1],4)*pow(-pow(m[0],2) + pow(m[1],2),5));
}

std::complex<double> LFh38(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],6) + 2*pow(m[0],4)*pow(m[1],2) + pow(m[0],2)*pow(m[1],4) - 2*pow(m[1],6) + 2*(pow(m[0],6) - 3*pow(m[0],4)*pow(m[1],2) + 3*pow(m[0],2)*pow(m[1],4))*log(mubarsq/pow(m[0],2)) - 2*pow(m[1],6)*log(mubarsq/pow(m[1],2)))/(2*pow(pow(m[0],2) - pow(m[1],2),3));
}

std::complex<double> LFh39(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*(pow(m[0],4) - 4*pow(m[0],2)*pow(m[1],2) + 3*pow(m[1],4) + 2*pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),3);
}

std::complex<double> LFh40(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],4) - pow(m[1],4) - 2*pow(m[0],2)*pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/(2*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),3));
}

std::complex<double> LFh41(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*(pow(m[0],6) - 6*pow(m[0],4)*pow(m[1],2) + 3*pow(m[0],2)*pow(m[1],4) + 2*pow(m[1],6) + 6*pow(m[0],2)*pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),4);
}

std::complex<double> LFh42(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],4) + 4*pow(m[0],2)*pow(m[1],2) - 5*pow(m[1],4) - 2*(2*pow(m[0],2)*pow(m[1],2) + pow(m[1],4))*log(pow(m[0],2)/pow(m[1],2)))/(2*pow(pow(m[0],2) - pow(m[1],2),4));
}

std::complex<double> LFh43(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*(-5*pow(m[0],4) + 4*pow(m[0],2)*pow(m[1],2) + pow(m[1],4) + 2*(pow(m[0],4) + 2*pow(m[0],2)*pow(m[1],2))*log(pow(m[0],2)/pow(m[1],2)))/(pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),4));
}

std::complex<double> LFh44(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],8) + 8*pow(m[0],6)*pow(m[1],2) - 8*pow(m[0],2)*pow(m[1],6) + pow(m[1],8) - 12*pow(m[0],4)*pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/(2*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh45(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],6) + 9*pow(m[0],4)*pow(m[1],2) - 9*pow(m[0],2)*pow(m[1],4) - pow(m[1],6) - 6*pow(m[0],2)*pow(m[1],2)*(pow(m[0],2) + pow(m[1],2))*log(pow(m[0],2)/pow(m[1],2)))/(2*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh46(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (3*(pow(m[0],4) - pow(m[1],4)) - (pow(m[0],4) + 4*pow(m[0],2)*pow(m[1],2) + pow(m[1],4))*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),5);
}

std::complex<double> LFh47(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],6) + 9*pow(m[0],4)*pow(m[1],2) - 9*pow(m[0],2)*pow(m[1],4) - pow(m[1],6) - 6*pow(m[0],2)*pow(m[1],2)*(pow(m[0],2) + pow(m[1],2))*log(pow(m[0],2)/pow(m[1],2)))/(2*pow(m[0],2)*pow(m[1],2)*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh48(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[0],6) + 9*pow(m[0],4)*pow(m[1],2) - 18*pow(m[0],2)*pow(m[1],4) + 11*pow(m[1],6) + 6*pow(m[1],6)*log(pow(m[0],2)/pow(m[1],2)))/(6*pow(pow(m[0],2) - pow(m[1],2),4));
}

std::complex<double> LFh49(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],6) - 6*pow(m[0],4)*pow(m[1],2) + 3*pow(m[0],2)*pow(m[1],4) + 2*pow(m[1],6) + 6*pow(m[0],2)*pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/(6*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),4));
}

std::complex<double> LFh50(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.16666666666666666*(2*pow(m[0],6) + 3*pow(m[0],4)*pow(m[1],2) - 6*pow(m[0],2)*pow(m[1],4) + pow(m[1],6) - 6*pow(m[0],4)*pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],3) - m[0]*pow(m[1],2),4);
}

std::complex<double> LFh51(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],8) + 6*pow(m[0],6)*pow(m[1],2) - 18*pow(m[0],4)*pow(m[1],4) + 10*pow(m[0],2)*pow(m[1],6) + 3*pow(m[1],8) + 12*pow(m[0],2)*pow(m[1],6)*log(pow(m[0],2)/pow(m[1],2)))/(3*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh52(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],6) - 9*pow(m[0],4)*pow(m[1],2) - 9*pow(m[0],2)*pow(m[1],4) + 17*pow(m[1],6) + 6*(3*pow(m[0],2)*pow(m[1],4) + pow(m[1],6))*log(pow(m[0],2)/pow(m[1],2)))/(6*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh53(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],6) - 9*pow(m[0],4)*pow(m[1],2) + 9*pow(m[0],2)*pow(m[1],4) + pow(m[1],6) + 6*pow(m[0],2)*pow(m[1],2)*(pow(m[0],2) + pow(m[1],2))*log(pow(m[0],2)/pow(m[1],2)))/(3*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh54(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-17*pow(m[0],6) + 9*pow(m[0],4)*pow(m[1],2) + 9*pow(m[0],2)*pow(m[1],4) - pow(m[1],6) + 6*(pow(m[0],6) + 3*pow(m[0],4)*pow(m[1],2))*log(pow(m[0],2)/pow(m[1],2)))/(6*pow(m[0],4)*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh55(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.08333333333333333*(3*pow(m[0],8) - 16*pow(m[0],6)*pow(m[1],2) + 36*pow(m[0],4)*pow(m[1],4) - 48*pow(m[0],2)*pow(m[1],6) + 25*pow(m[1],8) + 12*pow(m[1],8)*log(pow(m[0],2)/pow(m[1],2)))/pow(pow(m[0],2) - pow(m[1],2),5);
}

std::complex<double> LFh56(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],8) - 6*pow(m[0],6)*pow(m[1],2) + 18*pow(m[0],4)*pow(m[1],4) - 10*pow(m[0],2)*pow(m[1],6) - 3*pow(m[1],8) - 12*pow(m[0],2)*pow(m[1],6)*log(pow(m[0],2)/pow(m[1],2)))/(12*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh57(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-pow(m[0],8) + 8*pow(m[0],6)*pow(m[1],2) - 8*pow(m[0],2)*pow(m[1],6) + pow(m[1],8) - 12*pow(m[0],4)*pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/(12*pow(m[0],4)*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh58(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.08333333333333333*(-3*pow(m[0],8) - 10*pow(m[0],6)*pow(m[1],2) + 18*pow(m[0],4)*pow(m[1],4) - 6*pow(m[0],2)*pow(m[1],6) + pow(m[1],8) + 12*pow(m[0],6)*pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/(pow(m[0],6)*pow(pow(m[0],2) - pow(m[1],2),5));
}

std::complex<double> LFh59(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],4)*(pow(m[1],2) - pow(m[2],2))*log(mubarsq/pow(m[0],2)) + pow(m[1],4)*(-pow(m[0],2) + pow(m[2],2))*log(mubarsq/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2)) + pow(m[2],4)*log(mubarsq/pow(m[2],2))))/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh60(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],2)*(-pow(m[0],2) + pow(m[2],2))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2)))/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh61(const std::vector<std::complex<double> >& m, double mubarsq) {
    return ((pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/(-pow(m[0],2) + pow(m[1],2)) + (pow(m[2],2)*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2)) - (pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(-2*pow(m[1],2) + pow(m[2],2)))*log(pow(m[0],2)/pow(m[2],2))))/pow(pow(m[0],2) - pow(m[2],2),2))/pow(pow(m[1],2) - pow(m[2],2),2);
}

std::complex<double> LFh62(const std::vector<std::complex<double> >& m, double mubarsq) {
    return ((pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/(-pow(m[0],2) + pow(m[1],2)) + ((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2)) + (pow(m[0],2)*pow(m[1],2) - pow(m[2],4))*log(pow(m[0],2)/pow(m[2],2)))/pow(pow(m[0],2) - pow(m[2],2),2))/pow(pow(m[1],2) - pow(m[2],2),2);
}

std::complex<double> LFh63(const std::vector<std::complex<double> >& m, double mubarsq) {
    return ((2*pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/(-pow(m[0],2) + pow(m[1],2)) + (-((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[0],2)*(-3*pow(m[1],2) + pow(m[2],2)) + pow(m[2],2)*(pow(m[1],2) + pow(m[2],2)))) + 2*(pow(m[0],4)*pow(m[1],4) + pow(m[1],2)*pow(m[2],6) + pow(m[0],2)*(-3*pow(m[1],2)*pow(m[2],4) + pow(m[2],6)))*log(pow(m[0],2)/pow(m[2],2)))/pow(pow(m[0],2) - pow(m[2],2),3))/(2*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh64(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*pow(m[1],2)*pow(m[2],2)*pow(-pow(m[0],2) + pow(m[2],2),3)*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[2],2)*(pow(m[1],2) - 3*pow(m[2],2)) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))) - 2*pow(m[2],2)*(pow(m[0],4)*pow(m[1],2) + pow(m[2],6) + pow(m[0],2)*(pow(m[1],4) - 3*pow(m[1],2)*pow(m[2],2)))*log(pow(m[0],2)/pow(m[2],2))))/(2*(pow(m[0],2) - pow(m[1],2))*pow(m[2],2)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh65(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-6*pow(m[1],2)*pow(-(pow(m[0],2)*m[2]) + pow(m[2],3),4)*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(2*pow(m[1],4)*pow(m[2],4) - 7*pow(m[1],2)*pow(m[2],6) + 11*pow(m[2],8) + pow(m[0],4)*(-pow(m[1],4) + 5*pow(m[1],2)*pow(m[2],2) + 2*pow(m[2],4)) + pow(m[0],2)*(5*pow(m[1],4)*pow(m[2],2) - 10*pow(m[1],2)*pow(m[2],4) - 7*pow(m[2],6))) + 6*pow(m[2],4)*(pow(m[0],6)*pow(m[1],2) - pow(m[2],8) + pow(m[0],4)*(pow(m[1],4) - 4*pow(m[1],2)*pow(m[2],2)) + pow(m[0],2)*(pow(m[1],6) - 4*pow(m[1],4)*pow(m[2],2) + 6*pow(m[1],2)*pow(m[2],4)))*log(pow(m[0],2)/pow(m[2],2))))/(6*(pow(m[0],2) - pow(m[1],2))*pow(m[2],4)*pow(pow(m[0],2) - pow(m[2],2),4)*pow(pow(m[1],2) - pow(m[2],2),4));
}

std::complex<double> LFh66(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],2)*(-2*pow(m[0],2)*pow(m[2],4) + pow(m[1],2)*pow(m[2],4) - pow(m[0],4)*(pow(m[1],2) - 2*pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*(pow(m[1],2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2)) + (-pow(m[0],2) + pow(m[1],2))*pow(m[2],4)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh67(const std::vector<std::complex<double> >& m, double mubarsq) {
    return ((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],4) + pow(m[0],2)*pow(m[2],2))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2)) + (-pow(m[0],2) + pow(m[1],2))*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh68(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-(pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[1],4) - 2*pow(m[0],2)*pow(m[2],2) + pow(m[1],2)*pow(m[2],2))*log(pow(m[0],2)/pow(m[1],2))) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(-2*pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))) - (pow(m[0],2) - pow(m[1],2))*pow(m[2],2)*(2*pow(m[0],2)*pow(m[1],2) - pow(m[2],2)*(pow(m[1],2) + pow(m[2],2)))*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh69(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(pow(m[0],2) - pow(m[2],2),2)*(-2*pow(m[1],4) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(2*pow(m[0],2) - pow(m[1],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2)) + (pow(m[0],2) - pow(m[1],2))*(-2*pow(m[2],4) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2)))*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh70(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),3)) + 1/(2*(pow(m[0],2) - pow(m[2],2))*pow(-(pow(m[1],2)*m[2]) + pow(m[2],3),2)) + log(pow(m[0],2)/pow(m[1],2))/((-pow(m[0],2) + pow(m[1],2))*pow(pow(m[1],2) - pow(m[2],2),3)) - (pow(m[1],2)*(-3*pow(m[0],2) + 4*pow(m[1],2) - pow(m[2],2))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),4)) + ((2*pow(m[0],2) + pow(m[1],2) - 3*pow(m[2],2))*log(pow(m[0],2)/pow(m[2],2)))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(-pow(m[1],2) + pow(m[2],2),3)) - (pow(m[2],2)*(3*pow(m[0],4) + pow(m[1],4) - 4*pow(m[1],2)*pow(m[2],2) + 6*pow(m[2],4) + 2*pow(m[0],2)*(pow(m[1],2) - 4*pow(m[2],2)))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),4));
}

std::complex<double> LFh71(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],6)*pow(m[2],2) + pow(m[0],4)*pow(m[2],4) + pow(m[0],2)*(pow(m[1],6) - 3*pow(m[1],4)*pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[0],2)*(pow(m[1],2) - 3*pow(m[2],2)) + pow(m[1],2)*(pow(m[1],2) + pow(m[2],2))) + 2*pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],4)*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(pow(m[0],2) - pow(m[1],2),3)*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh72(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*pow(m[1],2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],6) - 3*pow(m[0],2)*pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*pow(m[2],2)*(pow(m[0],2) + pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(-3*pow(m[1],4) + pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))) - 2*pow(m[1],2)*pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(m[1],2)*pow(-pow(m[0],2) + pow(m[1],2),3)*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh73(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),2)) + 1/(2*(pow(m[0],2) - pow(m[1],2))*pow(pow(m[1],3) - m[1]*pow(m[2],2),2)) + ((2*pow(m[0],2) - 3*pow(m[1],2) + pow(m[2],2))*log(pow(m[0],2)/pow(m[1],2)))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),3)) + (pow(m[1],2)*(3*pow(m[0],4) + 6*pow(m[1],4) - 4*pow(m[1],2)*pow(m[2],2) + pow(m[2],4) + pow(m[0],2)*(-8*pow(m[1],2) + 2*pow(m[2],2)))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(-pow(m[0],2) + pow(m[1],2),3)*pow(pow(m[1],2) - pow(m[2],2),4)) + log(pow(m[0],2)/pow(m[2],2))/((-pow(m[0],2) + pow(m[2],2))*pow(-pow(m[1],2) + pow(m[2],2),3)) - (pow(m[2],2)*(-3*pow(m[0],2) - pow(m[1],2) + 4*pow(m[2],2))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),4));
}

std::complex<double> LFh74(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (6*pow(m[1],4)*(pow(m[0],2) - pow(m[2],2))*(-pow(m[1],8) + 6*pow(m[0],2)*pow(m[1],4)*pow(m[2],2) - 4*pow(m[0],2)*pow(m[1],2)*pow(m[2],2)*(pow(m[0],2) + pow(m[2],2)) + pow(m[0],2)*pow(m[2],2)*(pow(m[0],4) + pow(m[0],2)*pow(m[2],2) + pow(m[2],4)))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(-11*pow(m[1],8) + 7*pow(m[1],6)*pow(m[2],2) - 2*pow(m[1],4)*pow(m[2],4) + pow(m[0],4)*(-2*pow(m[1],4) - 5*pow(m[1],2)*pow(m[2],2) + pow(m[2],4)) + pow(m[0],2)*(7*pow(m[1],6) + 10*pow(m[1],4)*pow(m[2],2) - 5*pow(m[1],2)*pow(m[2],4))) + 6*pow(m[1],4)*pow(-pow(m[0],2) + pow(m[1],2),3)*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(6*pow(m[1],4)*pow(pow(m[0],2) - pow(m[1],2),4)*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[2],2),4));
}

std::complex<double> LFh75(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],4)*(pow(m[1],2) - pow(m[2],2))*(pow(m[0],4) + 3*pow(m[1],2)*pow(m[2],2) - 2*pow(m[0],2)*(pow(m[1],2) + pow(m[2],2)))*log(mubarsq/pow(m[0],2)) + pow(m[1],6)*pow(pow(m[0],2) - pow(m[2],2),2)*log(mubarsq/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(-(pow(m[1],2)*pow(m[2],2)) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))) + (pow(m[0],2) - pow(m[1],2))*pow(m[2],6)*log(mubarsq/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh76(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],4)*pow(pow(m[0],2) - pow(m[2],2),2)*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*(pow(m[0],2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2)) + (pow(m[0],2) - pow(m[1],2))*pow(m[2],4)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh77(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),2)*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2)) + (pow(m[0],2) - pow(m[1],2))*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh78(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],6)*pow(pow(m[0],2) - pow(m[2],2),3)*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[0],2)*pow(m[2],4) - pow(m[1],2)*pow(m[2],4) + pow(m[0],4)*(-pow(m[1],2) + pow(m[2],2))) + (pow(m[0],2) - pow(m[1],2))*pow(m[2],4)*(-(pow(m[1],2)*pow(m[2],2)) + pow(m[0],2)*(3*pow(m[1],2) - 2*pow(m[2],2)))*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh79(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],4)*pow(pow(m[0],2) - pow(m[2],2),3)*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(pow(m[1],2) - 2*pow(m[2],2))) + (pow(m[0],2) - pow(m[1],2))*pow(m[2],2)*(-pow(m[2],4) + pow(m[0],2)*(2*pow(m[1],2) - pow(m[2],2)))*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh80(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),3)*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - 2*pow(m[1],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2)) + (pow(m[0],2) - pow(m[1],2))*(pow(m[0],2)*pow(m[1],2) + pow(m[1],2)*pow(m[2],2) - 2*pow(m[2],4))*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh81(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)) + ((pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(-5*pow(m[1],4)*pow(m[2],2) + 9*pow(m[1],2)*pow(m[2],4) - 2*pow(m[2],6) + pow(m[0],4)*(pow(m[1],2) + pow(m[2],2)) - pow(m[0],2)*(pow(m[1],4) - 2*pow(m[1],2)*pow(m[2],2) + 5*pow(m[2],4))) - 2*(pow(m[0],2) - pow(m[1],2))*pow(m[2],2)*(pow(m[0],4)*pow(m[1],2) + 2*pow(m[0],2)*(pow(m[1],4) - 2*pow(m[1],2)*pow(m[2],2)) + pow(m[2],2)*(pow(m[1],4) - 3*pow(m[1],2)*pow(m[2],2) + 3*pow(m[2],4)))*log(pow(m[0],2)/pow(m[2],2))))/(pow(m[2],2)*pow(pow(m[0],2) - pow(m[2],2),4)))/(2*pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh82(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],4)*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(2*pow(m[1],2) - 3*pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(-(pow(m[0],2)*pow(m[1],4)) + pow(m[1],4)*pow(m[2],2) + pow(m[0],4)*(-pow(m[1],2) + pow(m[2],2))) - pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],6)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh83(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[1],4) + pow(m[0],2)*(pow(m[1],2) - 2*pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(-2*pow(m[1],2) + pow(m[2],2))) - pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],4)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh84(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -((pow(pow(m[0],2) - pow(m[2],2),2)*(-2*pow(m[1],4) + pow(m[0],2)*pow(m[2],2) + pow(m[1],2)*pow(m[2],2))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) + pow(m[1],2) - 2*pow(m[2],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2)) - pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)));
}

std::complex<double> LFh85(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),3)*(pow(m[1],4) - pow(m[0],2)*pow(m[2],2))*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[0],4)*(pow(m[1],2) + pow(m[2],2)) + pow(m[1],2)*pow(m[2],2)*(pow(m[1],2) + pow(m[2],2)) + pow(m[0],2)*(pow(m[1],4) - 6*pow(m[1],2)*pow(m[2],2) + pow(m[2],4))) - 2*pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],2)*(pow(m[0],2)*pow(m[1],2) - pow(m[2],4))*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh86(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-(pow(pow(m[0],2) - pow(m[2],2),3)*(-3*pow(m[1],4) + pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2))) + (pow(m[0],2) - pow(m[1],2))*(2*(pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[0],4) + pow(m[1],4) - pow(m[1],2)*pow(m[2],2) + pow(m[2],4) - pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))) + pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[2],2)*(pow(m[1],2) - 3*pow(m[2],2)) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2)))*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh87(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[1],8) + pow(m[0],4)*pow(m[2],4) + 2*pow(m[0],2)*pow(m[1],2)*(pow(m[1],4) - 3*pow(m[1],2)*pow(m[2],2) + pow(m[2],4)))*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(-3*pow(m[1],4)*pow(m[2],2) + pow(m[1],2)*pow(m[2],4) + pow(m[0],4)*(pow(m[1],2) - 3*pow(m[2],2)) + pow(m[0],2)*(5*pow(m[1],4) - 6*pow(m[1],2)*pow(m[2],2) + 5*pow(m[2],4))) + 2*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],4)*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh88(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*(3*pow(m[1],6) - 3*pow(m[1],4)*pow(m[2],2) + pow(m[0],2)*pow(m[2],2)*(pow(m[0],2) + 2*pow(m[2],2)) + pow(m[1],2)*(-4*pow(m[0],2)*pow(m[2],2) + pow(m[2],4)))*log(pow(m[0],2)/pow(m[1],2)) + ((pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(-2*pow(m[1],6) + 9*pow(m[1],4)*pow(m[2],2) - 5*pow(m[1],2)*pow(m[2],4) + pow(m[0],4)*(pow(m[1],2) + pow(m[2],2)) - pow(m[0],2)*(5*pow(m[1],4) - 2*pow(m[1],2)*pow(m[2],2) + pow(m[2],4))) + 2*pow(m[1],2)*pow(-pow(m[0],2) + pow(m[1],2),3)*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),2)))/(2*pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[1],2) - pow(m[2],2),3));
}

std::complex<double> LFh89(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*(2*pow(m[1],6)*pow(pow(m[0],2) - pow(m[2],2),3)*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*(pow(m[0],2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[0],4) + 5*pow(m[1],2)*pow(m[2],2) - 3*pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))) - 2*pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],6)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),3)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh90(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[1],4)*pow(pow(m[0],2) - pow(m[2],2),3)*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[0],4) - 3*pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))) + 2*pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],4)*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),3)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh91(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[0],2)*pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),3)*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(3*pow(m[0],4) - pow(m[1],2)*pow(m[2],2) - pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))) + 2*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),3)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh92(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[1],4)*pow(pow(m[0],2) - pow(m[2],2),4)*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(5*pow(m[1],4)*pow(m[2],2) - 3*pow(m[1],2)*pow(m[2],4) + pow(m[0],4)*(-3*pow(m[1],2) + 5*pow(m[2],2)) + pow(m[0],2)*(pow(m[1],4) - 6*pow(m[1],2)*pow(m[2],2) + pow(m[2],4))) + 2*pow(pow(m[0],2) - pow(m[1],2),2)*pow(m[2],2)*(pow(m[2],2)*(pow(m[1],2) - 2*pow(m[2],2)) + pow(m[0],2)*(2*pow(m[1],2) - pow(m[2],2)))*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),4)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh93(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[0],2)*pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),4)*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(2*pow(m[0],6) + pow(m[1],2)*pow(m[2],2)*(pow(m[1],2) - pow(m[2],2)) + pow(m[0],4)*(-9*pow(m[1],2) + 5*pow(m[2],2)) + pow(m[0],2)*(5*pow(m[1],4) - 2*pow(m[1],2)*pow(m[2],2) - pow(m[2],4))) + 2*pow(pow(m[0],3) - m[0]*pow(m[1],2),2)*(pow(m[0],2)*pow(m[1],2) + 2*pow(m[1],2)*pow(m[2],2) - 3*pow(m[2],4))*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),4)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh94(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*(2*pow(m[1],6)*pow(pow(m[0],2) - pow(m[2],2),3)*(pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(3*pow(m[1],2) - 4*pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(-4*pow(m[0],2)*pow(m[1],6)*pow(m[2],2) + 2*pow(m[1],6)*pow(m[2],4) + pow(m[0],8)*(-pow(m[1],2) + pow(m[2],2)) + pow(m[0],6)*(5*pow(m[1],4) - 2*pow(m[1],2)*pow(m[2],2) - 3*pow(m[2],4)) + pow(m[0],4)*(2*pow(m[1],6) - 7*pow(m[1],4)*pow(m[2],2) + 7*pow(m[1],2)*pow(m[2],4))) + 2*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],8)*log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh95(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[1],4)*pow(pow(m[0],2) - pow(m[2],2),3)*(pow(m[1],4) + pow(m[0],2)*(2*pow(m[1],2) - 3*pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(-2*pow(m[1],4)*pow(m[2],4) + pow(m[0],6)*(-pow(m[1],2) + pow(m[2],2)) + pow(m[0],4)*(-5*pow(m[1],4) + 2*pow(m[1],2)*pow(m[2],2) + pow(m[2],4)) + pow(m[0],2)*(9*pow(m[1],4)*pow(m[2],2) - 5*pow(m[1],2)*pow(m[2],4))) - 2*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],6)*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh96(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (-2*pow(m[1],2)*pow(pow(m[0],2) - pow(m[2],2),3)*(2*pow(m[1],4) - pow(m[1],2)*pow(m[2],2) + pow(m[0],2)*(pow(m[1],2) - 2*pow(m[2],2)))*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(-3*pow(m[1],4)*pow(m[2],2) + 5*pow(m[1],2)*pow(m[2],4) + pow(m[0],4)*(5*pow(m[1],2) - 3*pow(m[2],2)) + pow(m[0],2)*(pow(m[1],4) - 6*pow(m[1],2)*pow(m[2],2) + pow(m[2],4))) - 2*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],4)*log(pow(m[0],2)/pow(m[2],2))))/(2*pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh97(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*(-3*pow(m[1],4) + pow(m[0],2)*pow(m[2],2) + 2*pow(m[1],2)*pow(m[2],2))*log(pow(m[0],2)/pow(m[1],2)) + ((pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(2*pow(m[0],6) - pow(m[1],4)*pow(m[2],2) + pow(m[1],2)*pow(m[2],4) + pow(m[0],4)*(5*pow(m[1],2) - 9*pow(m[2],2)) - pow(m[0],2)*(pow(m[1],4) + 2*pow(m[1],2)*pow(m[2],2) - 5*pow(m[2],4))) - 2*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(pow(m[0],2)*pow(pow(m[0],2) - pow(m[2],2),3)))/(2*pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[1],2) - pow(m[2],2),2));
}

std::complex<double> LFh98(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (6*pow(m[1],8)*pow(pow(m[0],2) - pow(m[2],2),4)*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*(pow(m[0],2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(2*pow(m[0],8) + 26*pow(m[1],4)*pow(m[2],4) - 7*pow(m[0],6)*(pow(m[1],2) + pow(m[2],2)) - 31*pow(m[0],2)*pow(m[1],2)*pow(m[2],2)*(pow(m[1],2) + pow(m[2],2)) + pow(m[0],4)*(11*pow(m[1],4) + 26*pow(m[1],2)*pow(m[2],2) + 11*pow(m[2],4))) + 6*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],8)*log(pow(m[0],2)/pow(m[2],2))))/(6*pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[0],2) - pow(m[2],2),4)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh99(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (6*pow(m[1],6)*pow(pow(m[0],2) - pow(m[2],2),4)*log(pow(m[0],2)/pow(m[1],2)) + (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[0],8) - 11*pow(m[1],4)*pow(m[2],4) - 5*pow(m[0],6)*(pow(m[1],2) + pow(m[2],2)) + 7*pow(m[0],2)*pow(m[1],2)*pow(m[2],2)*(pow(m[1],2) + pow(m[2],2)) - 2*pow(m[0],4)*(pow(m[1],4) - 5*pow(m[1],2)*pow(m[2],2) + pow(m[2],4))) - 6*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],6)*log(pow(m[0],2)/pow(m[2],2))))/(6*pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[0],2) - pow(m[2],2),4)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh100(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (6*pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)) - ((pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(2*pow(m[0],8) + 2*pow(m[1],4)*pow(m[2],4) + 5*pow(m[0],6)*(pow(m[1],2) + pow(m[2],2)) + 5*pow(m[0],2)*pow(m[1],2)*pow(m[2],2)*(pow(m[1],2) + pow(m[2],2)) - pow(m[0],4)*(pow(m[1],4) + 22*pow(m[1],2)*pow(m[2],2) + pow(m[2],4))) + 6*pow(m[0],2)*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],4)*log(pow(m[0],2)/pow(m[2],2))))/(pow(m[0],2)*pow(pow(m[0],2) - pow(m[2],2),4)))/(6*pow(pow(m[0],2) - pow(m[1],2),4)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh101(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (6*pow(m[1],2)*pow(pow(m[0],3) - m[0]*pow(m[2],2),4)*log(pow(m[0],2)/pow(m[1],2)) - (pow(m[0],2) - pow(m[1],2))*((pow(m[0],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[2],2))*(11*pow(m[0],8) - pow(m[1],4)*pow(m[2],4) - 7*pow(m[0],6)*(pow(m[1],2) + pow(m[2],2)) + 5*pow(m[0],2)*pow(m[1],2)*pow(m[2],2)*(pow(m[1],2) + pow(m[2],2)) + 2*pow(m[0],4)*(pow(m[1],4) - 5*pow(m[1],2)*pow(m[2],2) + pow(m[2],4))) + 6*pow(m[0],4)*pow(pow(m[0],2) - pow(m[1],2),3)*pow(m[2],2)*log(pow(m[0],2)/pow(m[2],2))))/(6*pow(m[0],4)*pow(pow(m[0],2) - pow(m[1],2),4)*pow(pow(m[0],2) - pow(m[2],2),4)*(pow(m[1],2) - pow(m[2],2)));
}

std::complex<double> LFh102(const std::vector<std::complex<double> >& m, double mubarsq) {
    return ((pow(m[1],4)*log(pow(m[0],2)/pow(m[1],2)))/(-pow(m[0],2) + pow(m[1],2)) + ((pow(m[2],4)*(pow(m[1],2) - pow(m[3],2))*log(pow(m[0],2)/pow(m[2],2)))/(pow(m[0],2) - pow(m[2],2)) + ((pow(m[1],2) - pow(m[2],2))*pow(m[3],4)*log(pow(m[0],2)/pow(m[3],2)))/(-pow(m[0],2) + pow(m[3],2)))/(pow(m[2],2) - pow(m[3],2)))/((pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2)));
}

std::complex<double> LFh103(const std::vector<std::complex<double> >& m, double mubarsq) {
    return ((pow(m[1],2)*log(pow(m[0],2)/pow(m[1],2)))/(-pow(m[0],2) + pow(m[1],2)) + ((pow(m[2],2)*(pow(m[1],2) - pow(m[3],2))*log(pow(m[0],2)/pow(m[2],2)))/(pow(m[0],2) - pow(m[2],2)) + ((pow(m[1],2) - pow(m[2],2))*pow(m[3],2)*log(pow(m[0],2)/pow(m[3],2)))/(-pow(m[0],2) + pow(m[3],2)))/(pow(m[2],2) - pow(m[3],2)))/((pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2)));
}

std::complex<double> LFh104(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],4)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[0],2) - pow(m[3],2),2)) + (pow(m[1],4)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[3],2),2)) + (pow(m[2],4)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*pow(pow(m[2],2) - pow(m[3],2),2)) - (pow(m[3],4)*(-2*pow(m[2],2)*pow(m[3],2) + 3*pow(m[3],4) + pow(m[1],2)*(pow(m[2],2) - 2*pow(m[3],2)) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2) - 2*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)) + (pow(m[3],2)*(1 + 2*log(pow(m[0],2)/pow(m[3],2))))/((-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh105(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[0],2) - pow(m[3],2),2)) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[3],2),2)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*pow(pow(m[2],2) - pow(m[3],2),2)) + log(pow(m[0],2)/pow(m[3],2))/((-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))) - (pow(m[3],2)*(-2*pow(m[2],2)*pow(m[3],2) + 3*pow(m[3],4) + pow(m[1],2)*(pow(m[2],2) - 2*pow(m[3],2)) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2) - 2*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[3],2),2));
}

std::complex<double> LFh106(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[0],2) - pow(m[3],2),3)) - 1/(2*pow(m[3],2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[3],2),3)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*pow(pow(m[2],2) - pow(m[3],2),3)) - ((-2*pow(m[2],2)*pow(m[3],2) + 3*pow(m[3],4) + pow(m[1],2)*(pow(m[2],2) - 2*pow(m[3],2)) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2) - 2*pow(m[3],2)))*log(pow(m[0],2)/pow(m[3],2)))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)) - (pow(m[3],2)*(3*pow(m[2],4)*pow(m[3],4) - 8*pow(m[2],2)*pow(m[3],6) + 6*pow(m[3],8) + pow(m[1],4)*(pow(m[2],4) - 3*pow(m[2],2)*pow(m[3],2) + 3*pow(m[3],4)) + pow(m[1],2)*(-3*pow(m[2],4)*pow(m[3],2) + 9*pow(m[2],2)*pow(m[3],4) - 8*pow(m[3],6)) + pow(m[0],4)*(pow(m[1],4) + pow(m[2],4) - 3*pow(m[2],2)*pow(m[3],2) + 3*pow(m[3],4) + pow(m[1],2)*(pow(m[2],2) - 3*pow(m[3],2))) + pow(m[0],2)*(-3*pow(m[2],4)*pow(m[3],2) + 9*pow(m[2],2)*pow(m[3],4) - 8*pow(m[3],6) + pow(m[1],4)*(pow(m[2],2) - 3*pow(m[3],2)) + pow(m[1],2)*pow(pow(m[2],2) - 3*pow(m[3],2),2)))*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),3)*pow(-pow(m[1],2) + pow(m[3],2),3)*pow(-pow(m[2],2) + pow(m[3],2),3));
}

std::complex<double> LFh107(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],4)/((pow(m[0],2) - pow(m[1],2))*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[0],2) - pow(m[3],2))) + (pow(m[1],4)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*pow(pow(m[1],2) - pow(m[2],2),2)*(pow(m[1],2) - pow(m[3],2))) - (pow(m[2],4)*(3*pow(m[2],4) - 2*pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(-2*pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(pow(m[1],2) - 2*pow(m[2],2) + pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)) + (pow(m[2],2)*(1 + 2*log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],4)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[2],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2)));
}

std::complex<double> LFh108(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[0],2) - pow(m[3],2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*pow(pow(m[1],2) - pow(m[2],2),2)*(pow(m[1],2) - pow(m[3],2))) + log(pow(m[0],2)/pow(m[2],2))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) - (pow(m[2],2)*(3*pow(m[2],4) - 2*pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(-2*pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(pow(m[1],2) - 2*pow(m[2],2) + pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[2],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2)));
}

std::complex<double> LFh109(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),2)) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)) + log(pow(m[0],2)/pow(m[2],2))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*pow(pow(m[2],2) - pow(m[3],2),2)) - (pow(m[2],2)*(4*pow(m[2],4) - 2*pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(-3*pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(2*pow(m[1],2) - 3*pow(m[2],2) + pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[2],2) - pow(m[3],2),3)) + log(pow(m[0],2)/pow(m[3],2))/(pow(pow(m[2],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))) - (pow(m[3],2)*(-2*pow(m[2],2)*pow(m[3],2) + 4*pow(m[3],4) + pow(m[1],2)*(pow(m[2],2) - 3*pow(m[3],2)) + pow(m[0],2)*(2*pow(m[1],2) + pow(m[2],2) - 3*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)*pow(-pow(m[2],2) + pow(m[3],2),3));
}

std::complex<double> LFh110(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*pow(pow(m[0],2) - pow(m[2],2),3)*(pow(m[0],2) - pow(m[3],2))) - 1/(2*pow(m[2],2)*(-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*pow(pow(m[1],2) - pow(m[2],2),3)*(pow(m[1],2) - pow(m[3],2))) - ((3*pow(m[2],4) - 2*pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(-2*pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(pow(m[1],2) - 2*pow(m[2],2) + pow(m[3],2)))*log(pow(m[0],2)/pow(m[2],2)))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)) - (pow(m[2],2)*(6*pow(m[2],8) - 8*pow(m[2],6)*pow(m[3],2) + 3*pow(m[2],4)*pow(m[3],4) + pow(m[1],4)*(3*pow(m[2],4) - 3*pow(m[2],2)*pow(m[3],2) + pow(m[3],4)) + pow(m[1],2)*(-8*pow(m[2],6) + 9*pow(m[2],4)*pow(m[3],2) - 3*pow(m[2],2)*pow(m[3],4)) + pow(m[0],4)*(pow(m[1],4) + 3*pow(m[2],4) - 3*pow(m[2],2)*pow(m[3],2) + pow(m[3],4) + pow(m[1],2)*(-3*pow(m[2],2) + pow(m[3],2))) + pow(m[0],2)*(-8*pow(m[2],6) + 9*pow(m[2],4)*pow(m[3],2) - 3*pow(m[2],2)*pow(m[3],4) + pow(m[1],4)*(-3*pow(m[2],2) + pow(m[3],2)) + pow(m[1],2)*pow(-3*pow(m[2],2) + pow(m[3],2),2)))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),3)*pow(-pow(m[1],2) + pow(m[2],2),3)*pow(pow(m[2],2) - pow(m[3],2),3)) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/((-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*pow(-pow(m[2],2) + pow(m[3],2),3));
}

std::complex<double> LFh111(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],4)/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))) - (pow(m[1],4)*(3*pow(m[1],4) + pow(m[2],2)*pow(m[3],2) - 2*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(-2*pow(m[1],2) + pow(m[2],2) + pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)) + (pow(m[1],2)*(1 + 2*log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + (pow(m[2],4)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[1],2) - pow(m[2],2),2)*(-pow(m[0],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],4)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[1],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh112(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))) + log(pow(m[0],2)/pow(m[1],2))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) - (pow(m[1],2)*(3*pow(m[1],4) + pow(m[2],2)*pow(m[3],2) - 2*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(-2*pow(m[1],2) + pow(m[2],2) + pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[1],2) - pow(m[2],2),2)*(-pow(m[0],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[1],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh113(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[0],2) - pow(m[3],2),2)) + log(pow(m[0],2)/pow(m[1],2))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[3],2),2)) - (pow(m[1],2)*(4*pow(m[1],4) + pow(m[2],2)*pow(m[3],2) + pow(m[0],2)*(-3*pow(m[1],2) + 2*pow(m[2],2) + pow(m[3],2)) - pow(m[1],2)*(3*pow(m[2],2) + 2*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[3],2),3)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[1],2) - pow(m[2],2),2)*(-pow(m[0],2) + pow(m[2],2))*pow(pow(m[2],2) - pow(m[3],2),2)) + log(pow(m[0],2)/pow(m[3],2))/(pow(pow(m[1],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))) - (pow(m[3],2)*(-3*pow(m[2],2)*pow(m[3],2) + 4*pow(m[3],4) + pow(m[0],2)*(pow(m[1],2) + 2*pow(m[2],2) - 3*pow(m[3],2)) + pow(m[1],2)*(pow(m[2],2) - 2*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)*pow(-pow(m[1],2) + pow(m[3],2),3));
}

std::complex<double> LFh114(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[0],2) - pow(m[3],2))) + log(pow(m[0],2)/pow(m[1],2))/((-pow(m[0],2) + pow(m[1],2))*pow(pow(m[1],2) - pow(m[2],2),2)*(pow(m[1],2) - pow(m[3],2))) - (pow(m[1],2)*(4*pow(m[1],4) + pow(m[2],2)*pow(m[3],2) + pow(m[0],2)*(-3*pow(m[1],2) + pow(m[2],2) + 2*pow(m[3],2)) - pow(m[1],2)*(2*pow(m[2],2) + 3*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[3],2),2)) + log(pow(m[0],2)/pow(m[2],2))/(pow(pow(m[1],2) - pow(m[2],2),2)*(-pow(m[0],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) - (pow(m[2],2)*(4*pow(m[2],4) - 3*pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(-2*pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(pow(m[1],2) - 3*pow(m[2],2) + 2*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(-pow(m[1],2) + pow(m[2],2),3)*pow(pow(m[2],2) - pow(m[3],2),2)) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[1],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2)));
}

std::complex<double> LFh115(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/(pow(pow(m[0],2) - pow(m[1],2),3)*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))) - 1/(2*pow(m[1],2)*(-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + ((-3*pow(m[1],4) - pow(m[2],2)*pow(m[3],2) + pow(m[0],2)*(2*pow(m[1],2) - pow(m[2],2) - pow(m[3],2)) + 2*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)))*log(pow(m[0],2)/pow(m[1],2)))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)) - (pow(m[1],2)*(6*pow(m[1],8) + pow(m[2],4)*pow(m[3],4) - 8*pow(m[1],6)*(pow(m[2],2) + pow(m[3],2)) - 3*pow(m[1],2)*pow(m[2],2)*pow(m[3],2)*(pow(m[2],2) + pow(m[3],2)) + 3*pow(m[1],4)*(pow(m[2],4) + 3*pow(m[2],2)*pow(m[3],2) + pow(m[3],4)) + pow(m[0],4)*(3*pow(m[1],4) + pow(m[2],4) + pow(m[2],2)*pow(m[3],2) + pow(m[3],4) - 3*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2))) + pow(m[0],2)*(-8*pow(m[1],6) + 9*pow(m[1],4)*(pow(m[2],2) + pow(m[3],2)) + pow(m[2],2)*pow(m[3],2)*(pow(m[2],2) + pow(m[3],2)) - 3*pow(m[1],2)*pow(pow(m[2],2) + pow(m[3],2),2)))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[1],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[3],2),3)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*pow(-pow(m[1],2) + pow(m[2],2),3)*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/((-pow(m[0],2) + pow(m[3],2))*pow(-pow(m[1],2) + pow(m[3],2),3)*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh116(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (2*pow(m[0],4))/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))) - (pow(m[0],6)*(3*pow(m[0],4) + pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) - 2*pow(m[0],2)*(pow(m[1],2) + pow(m[2],2) + pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),2)) + (pow(m[1],6)*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + (pow(m[2],6)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],6)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh117(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))) - (pow(m[0],4)*(3*pow(m[0],4) + pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) - 2*pow(m[0],2)*(pow(m[1],2) + pow(m[2],2) + pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),2)) + (pow(m[1],4)*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + (pow(m[2],4)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],4)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh118(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -((pow(m[0],2)*(3*pow(m[0],4) + pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) - 2*pow(m[0],2)*(pow(m[1],2) + pow(m[2],2) + pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh119(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -((pow(m[0],2)*(4*pow(m[0],4) + pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(2*pow(m[2],2) + pow(m[3],2)) - pow(m[0],2)*(3*pow(m[1],2) + 3*pow(m[2],2) + 2*pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),3))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[1],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[3],2),2)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*(-pow(m[1],2) + pow(m[2],2))*pow(pow(m[2],2) - pow(m[3],2),2)) + log(pow(m[0],2)/pow(m[3],2))/(pow(pow(m[0],2) - pow(m[3],2),2)*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))) + (pow(m[3],2)*(-3*pow(m[2],2)*pow(m[3],2) + 4*pow(m[3],4) + pow(m[1],2)*(2*pow(m[2],2) - 3*pow(m[3],2)) + pow(m[0],2)*(pow(m[1],2) + pow(m[2],2) - 2*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),3)*pow(pow(m[1],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[3],2),2));
}

std::complex<double> LFh120(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -((pow(m[0],2)*(4*pow(m[0],4) + pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(pow(m[2],2) + 2*pow(m[3],2)) - pow(m[0],2)*(3*pow(m[1],2) + 2*pow(m[2],2) + 3*pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[0],2) - pow(m[3],2),2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*(pow(m[1],2) - pow(m[3],2))) + log(pow(m[0],2)/pow(m[2],2))/(pow(pow(m[0],2) - pow(m[2],2),2)*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[2],2)*(4*pow(m[2],4) - 3*pow(m[2],2)*pow(m[3],2) + pow(m[0],2)*(pow(m[1],2) - 2*pow(m[2],2) + pow(m[3],2)) + pow(m[1],2)*(-3*pow(m[2],2) + 2*pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)*(-pow(m[1],2) + pow(m[3],2)));
}

std::complex<double> LFh121(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))) - (pow(m[0],4)*(4*pow(m[0],4) + 2*pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) - pow(m[0],2)*(2*pow(m[1],2) + 3*(pow(m[2],2) + pow(m[3],2)))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),2)) + (pow(m[1],4)*(4*pow(m[1],4) + 2*pow(m[2],2)*pow(m[3],2) - 3*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(-2*pow(m[1],2) + pow(m[2],2) + pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)) + (pow(m[1],2)*(1 + 2*log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + (pow(m[2],4)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],4)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh122(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -((pow(m[0],2)*(4*pow(m[0],4) + 2*pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) - pow(m[0],2)*(2*pow(m[1],2) + 3*(pow(m[2],2) + pow(m[3],2)))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),2))) + log(pow(m[0],2)/pow(m[1],2))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + (pow(m[1],2)*(4*pow(m[1],4) + 2*pow(m[2],2)*pow(m[3],2) - 3*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) + pow(m[0],2)*(-2*pow(m[1],2) + pow(m[2],2) + pow(m[3],2)))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh123(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*1/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))) - (pow(m[0],2)*(3*pow(m[0],4) + pow(m[2],2)*pow(m[3],2) + pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)) - 2*pow(m[0],2)*(pow(m[1],2) + pow(m[2],2) + pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),2)) + (pow(m[0],4)*(6*pow(m[0],8) + pow(m[2],4)*pow(m[3],4) + pow(m[1],2)*pow(m[2],2)*pow(m[3],2)*(pow(m[2],2) + pow(m[3],2)) - 3*pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))*(pow(m[1],2) + pow(m[3],2))*(pow(m[2],2) + pow(m[3],2)) - 8*pow(m[0],6)*(pow(m[1],2) + pow(m[2],2) + pow(m[3],2)) + pow(m[1],4)*(pow(m[2],4) + pow(m[2],2)*pow(m[3],2) + pow(m[3],4)) + 3*pow(m[0],4)*(pow(m[1],4) + pow(m[2],4) + 3*pow(m[2],2)*pow(m[3],2) + pow(m[3],4) + 3*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[0],2) - pow(m[3],2),3)) + (pow(m[1],4)*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(-pow(m[0],2) + pow(m[1],2),3)*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + (pow(m[2],4)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(-pow(m[0],2) + pow(m[2],2),3)*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],4)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(-pow(m[0],2) + pow(m[3],2),3)*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh124(const std::vector<std::complex<double> >& m, double mubarsq) {
    return -0.5*1/(pow(m[0],2)*(pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))) + (pow(m[0],2)*(6*pow(m[0],8) + pow(m[2],4)*pow(m[3],4) + pow(m[1],2)*pow(m[2],2)*pow(m[3],2)*(pow(m[2],2) + pow(m[3],2)) - 3*pow(m[0],2)*(pow(m[1],2) + pow(m[2],2))*(pow(m[1],2) + pow(m[3],2))*(pow(m[2],2) + pow(m[3],2)) - 8*pow(m[0],6)*(pow(m[1],2) + pow(m[2],2) + pow(m[3],2)) + pow(m[1],4)*(pow(m[2],4) + pow(m[2],2)*pow(m[3],2) + pow(m[3],4)) + 3*pow(m[0],4)*(pow(m[1],4) + pow(m[2],4) + 3*pow(m[2],2)*pow(m[3],2) + pow(m[3],4) + 3*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2)))))/(pow(pow(m[0],2) - pow(m[1],2),3)*pow(pow(m[0],2) - pow(m[2],2),3)*pow(pow(m[0],2) - pow(m[3],2),3)) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(-pow(m[0],2) + pow(m[1],2),3)*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(-pow(m[0],2) + pow(m[2],2),3)*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(-pow(m[0],2) + pow(m[3],2),3)*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2)));
}

std::complex<double> LFh125(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],4)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))*(pow(m[0],2) - pow(m[4],2))) + (pow(m[1],4)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))*(pow(m[1],2) - pow(m[4],2))) + (pow(m[2],4)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))*(pow(m[2],2) - pow(m[4],2))) + (pow(m[3],4)*(1 + log(pow(m[0],2)/pow(m[3],2))))/((-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))*(pow(m[3],2) - pow(m[4],2))) + (pow(m[4],4)*(1 + log(pow(m[0],2)/pow(m[4],2))))/((-pow(m[0],2) + pow(m[4],2))*(-pow(m[1],2) + pow(m[4],2))*(-pow(m[2],2) + pow(m[4],2))*(-pow(m[3],2) + pow(m[4],2)));
}

std::complex<double> LFh126(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))*(pow(m[0],2) - pow(m[4],2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))*(pow(m[1],2) - pow(m[4],2))) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))*(pow(m[2],2) - pow(m[4],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/((-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))*(pow(m[3],2) - pow(m[4],2))) + (pow(m[4],2)*(1 + log(pow(m[0],2)/pow(m[4],2))))/((-pow(m[0],2) + pow(m[4],2))*(-pow(m[1],2) + pow(m[4],2))*(-pow(m[2],2) + pow(m[4],2))*(-pow(m[3],2) + pow(m[4],2)));
}

std::complex<double> LFh127(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))*pow(pow(m[0],2) - pow(m[4],2),2)) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))*pow(pow(m[1],2) - pow(m[4],2),2)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))*pow(pow(m[2],2) - pow(m[4],2),2)) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/((-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))*pow(pow(m[3],2) - pow(m[4],2),2)) + log(pow(m[0],2)/pow(m[4],2))/((-pow(m[0],2) + pow(m[4],2))*(-pow(m[1],2) + pow(m[4],2))*(-pow(m[2],2) + pow(m[4],2))*(-pow(m[3],2) + pow(m[4],2))) - (pow(m[4],2)*(2*pow(m[2],2)*pow(m[3],2)*pow(m[4],2) - 3*pow(m[2],2)*pow(m[4],4) - 3*pow(m[3],2)*pow(m[4],4) + 4*pow(m[4],6) + pow(m[1],2)*(2*pow(m[3],2)*pow(m[4],2) - 3*pow(m[4],4) - pow(m[2],2)*(pow(m[3],2) - 2*pow(m[4],2))) - pow(m[0],2)*(-2*pow(m[3],2)*pow(m[4],2) + 3*pow(m[4],4) + pow(m[2],2)*(pow(m[3],2) - 2*pow(m[4],2)) + pow(m[1],2)*(pow(m[2],2) + pow(m[3],2) - 2*pow(m[4],2))))*(1 + log(pow(m[0],2)/pow(m[4],2))))/(pow(pow(m[0],2) - pow(m[4],2),2)*pow(pow(m[1],2) - pow(m[4],2),2)*pow(pow(m[2],2) - pow(m[4],2),2)*pow(pow(m[3],2) - pow(m[4],2),2));
}

std::complex<double> LFh128(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*pow(pow(m[0],2) - pow(m[3],2),2)*(pow(m[0],2) - pow(m[4],2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*pow(pow(m[1],2) - pow(m[3],2),2)*(pow(m[1],2) - pow(m[4],2))) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*pow(pow(m[2],2) - pow(m[3],2),2)*(pow(m[2],2) - pow(m[4],2))) + log(pow(m[0],2)/pow(m[3],2))/((-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))*(pow(m[3],2) - pow(m[4],2))) - (pow(m[3],2)*(-3*pow(m[2],2)*pow(m[3],4) + 4*pow(m[3],6) + 2*pow(m[2],2)*pow(m[3],2)*pow(m[4],2) - 3*pow(m[3],4)*pow(m[4],2) + pow(m[1],2)*(-3*pow(m[3],4) + 2*pow(m[3],2)*pow(m[4],2) + pow(m[2],2)*(2*pow(m[3],2) - pow(m[4],2))) - pow(m[0],2)*(3*pow(m[3],4) - 2*pow(m[3],2)*pow(m[4],2) + pow(m[2],2)*(-2*pow(m[3],2) + pow(m[4],2)) + pow(m[1],2)*(pow(m[2],2) - 2*pow(m[3],2) + pow(m[4],2))))*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)*pow(pow(m[3],2) - pow(m[4],2),2)) + (pow(m[4],2)*(1 + log(pow(m[0],2)/pow(m[4],2))))/(pow(pow(m[3],2) - pow(m[4],2),2)*(-pow(m[0],2) + pow(m[4],2))*(-pow(m[1],2) + pow(m[4],2))*(-pow(m[2],2) + pow(m[4],2)));
}

std::complex<double> LFh129(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*pow(pow(m[0],2) - pow(m[2],2),2)*(pow(m[0],2) - pow(m[3],2))*(pow(m[0],2) - pow(m[4],2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*pow(pow(m[1],2) - pow(m[2],2),2)*(pow(m[1],2) - pow(m[3],2))*(pow(m[1],2) - pow(m[4],2))) + log(pow(m[0],2)/pow(m[2],2))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))*(pow(m[2],2) - pow(m[4],2))) - (pow(m[2],2)*(4*pow(m[2],6) - 3*pow(m[2],4)*pow(m[3],2) - 3*pow(m[2],4)*pow(m[4],2) + 2*pow(m[2],2)*pow(m[3],2)*pow(m[4],2) + pow(m[1],2)*(-3*pow(m[2],4) - pow(m[3],2)*pow(m[4],2) + 2*pow(m[2],2)*(pow(m[3],2) + pow(m[4],2))) + pow(m[0],2)*(-3*pow(m[2],4) - pow(m[3],2)*pow(m[4],2) + pow(m[1],2)*(2*pow(m[2],2) - pow(m[3],2) - pow(m[4],2)) + 2*pow(m[2],2)*(pow(m[3],2) + pow(m[4],2))))*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[2],2) - pow(m[3],2),2)*pow(pow(m[2],2) - pow(m[4],2),2)) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[2],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(pow(m[3],2) - pow(m[4],2))) + (pow(m[4],2)*(1 + log(pow(m[0],2)/pow(m[4],2))))/(pow(pow(m[2],2) - pow(m[4],2),2)*(-pow(m[0],2) + pow(m[4],2))*(-pow(m[1],2) + pow(m[4],2))*(-pow(m[3],2) + pow(m[4],2)));
}

std::complex<double> LFh130(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))*(pow(m[0],2) - pow(m[4],2))) + log(pow(m[0],2)/pow(m[1],2))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))*(pow(m[1],2) - pow(m[4],2))) - (pow(m[1],2)*(4*pow(m[1],6) - pow(m[2],2)*pow(m[3],2)*pow(m[4],2) - 3*pow(m[1],4)*(pow(m[2],2) + pow(m[3],2) + pow(m[4],2)) + 2*pow(m[1],2)*(pow(m[3],2)*pow(m[4],2) + pow(m[2],2)*(pow(m[3],2) + pow(m[4],2))) - pow(m[0],2)*(3*pow(m[1],4) + pow(m[3],2)*pow(m[4],2) + pow(m[2],2)*(pow(m[3],2) + pow(m[4],2)) - 2*pow(m[1],2)*(pow(m[2],2) + pow(m[3],2) + pow(m[4],2))))*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[1],2) - pow(m[2],2),2)*pow(pow(m[1],2) - pow(m[3],2),2)*pow(pow(m[1],2) - pow(m[4],2),2)) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[1],2) - pow(m[2],2),2)*(-pow(m[0],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))*(pow(m[2],2) - pow(m[4],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[1],2) - pow(m[3],2),2)*(-pow(m[0],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))*(pow(m[3],2) - pow(m[4],2))) + (pow(m[4],2)*(1 + log(pow(m[0],2)/pow(m[4],2))))/(pow(pow(m[1],2) - pow(m[4],2),2)*(-pow(m[0],2) + pow(m[4],2))*(-pow(m[2],2) + pow(m[4],2))*(-pow(m[3],2) + pow(m[4],2)));
}

std::complex<double> LFh131(const std::vector<std::complex<double> >& m, double mubarsq) {
    return (pow(m[0],2)*(-4*pow(m[0],6) + pow(m[2],2)*pow(m[3],2)*pow(m[4],2) + 3*pow(m[0],4)*(pow(m[1],2) + pow(m[2],2) + pow(m[3],2) + pow(m[4],2)) + pow(m[1],2)*(pow(m[3],2)*pow(m[4],2) + pow(m[2],2)*(pow(m[3],2) + pow(m[4],2))) - 2*pow(m[0],2)*(pow(m[3],2)*pow(m[4],2) + pow(m[2],2)*(pow(m[3],2) + pow(m[4],2)) + pow(m[1],2)*(pow(m[2],2) + pow(m[3],2) + pow(m[4],2)))))/(pow(pow(m[0],2) - pow(m[1],2),2)*pow(pow(m[0],2) - pow(m[2],2),2)*pow(pow(m[0],2) - pow(m[3],2),2)*pow(pow(m[0],2) - pow(m[4],2),2)) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/(pow(pow(m[0],2) - pow(m[1],2),2)*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))*(pow(m[1],2) - pow(m[4],2))) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/(pow(pow(m[0],2) - pow(m[2],2),2)*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))*(pow(m[2],2) - pow(m[4],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/(pow(pow(m[0],2) - pow(m[3],2),2)*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))*(pow(m[3],2) - pow(m[4],2))) + (pow(m[4],2)*(1 + log(pow(m[0],2)/pow(m[4],2))))/(pow(pow(m[0],2) - pow(m[4],2),2)*(-pow(m[1],2) + pow(m[4],2))*(-pow(m[2],2) + pow(m[4],2))*(-pow(m[3],2) + pow(m[4],2)));
}

std::complex<double> LFh132(const std::vector<std::complex<double> >& m, double mubarsq) {
    return pow(m[0],2)/((pow(m[0],2) - pow(m[1],2))*(pow(m[0],2) - pow(m[2],2))*(pow(m[0],2) - pow(m[3],2))*(pow(m[0],2) - pow(m[4],2))*(pow(m[0],2) - pow(m[5],2))) + (pow(m[1],2)*(1 + log(pow(m[0],2)/pow(m[1],2))))/((-pow(m[0],2) + pow(m[1],2))*(pow(m[1],2) - pow(m[2],2))*(pow(m[1],2) - pow(m[3],2))*(pow(m[1],2) - pow(m[4],2))*(pow(m[1],2) - pow(m[5],2))) + (pow(m[2],2)*(1 + log(pow(m[0],2)/pow(m[2],2))))/((-pow(m[0],2) + pow(m[2],2))*(-pow(m[1],2) + pow(m[2],2))*(pow(m[2],2) - pow(m[3],2))*(pow(m[2],2) - pow(m[4],2))*(pow(m[2],2) - pow(m[5],2))) + (pow(m[3],2)*(1 + log(pow(m[0],2)/pow(m[3],2))))/((-pow(m[0],2) + pow(m[3],2))*(-pow(m[1],2) + pow(m[3],2))*(-pow(m[2],2) + pow(m[3],2))*(pow(m[3],2) - pow(m[4],2))*(pow(m[3],2) - pow(m[5],2))) + (pow(m[4],2)*(1 + log(pow(m[0],2)/pow(m[4],2))))/((-pow(m[0],2) + pow(m[4],2))*(-pow(m[1],2) + pow(m[4],2))*(-pow(m[2],2) + pow(m[4],2))*(-pow(m[3],2) + pow(m[4],2))*(pow(m[4],2) - pow(m[5],2))) + (pow(m[5],2)*(1 + log(pow(m[0],2)/pow(m[5],2))))/((-pow(m[0],2) + pow(m[5],2))*(-pow(m[1],2) + pow(m[5],2))*(-pow(m[2],2) + pow(m[5],2))*(-pow(m[3],2) + pow(m[5],2))*(-pow(m[4],2) + pow(m[5],2)));
}


