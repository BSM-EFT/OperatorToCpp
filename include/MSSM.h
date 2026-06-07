#pragma once
#include <vector>
#include <string>
#include <complex>
#include <unordered_map>
#include <map>
#include <functional>
#include <cmath>

const double pi = 3.14159265;

struct Task {
    std::string name;
    std::function<std::complex<double>()> work;
};

class MSSM {
    private:
        double hbar = 1/(16 * pow(pi,2));
        double mubarsq = 1.0;
        std::complex<double> cgamma = 0.0;
        std::complex<double> g1 = 0.0;
        std::complex<double> g2 = 0.0;
        std::complex<double> g3 = 0.0;
        std::complex<double> m1 = 0.0;
        std::complex<double> m2 = 0.0;
        std::complex<double> m3 = 0.0;
        std::complex<double> mPhi = 0.0;
        std::complex<double> muTilde = 0.0;

        std::vector<std::complex<double> > mdt = {0.0, 0.0, 0.0};
        std::vector<std::complex<double> > met = {0.0, 0.0, 0.0};
        std::vector<std::complex<double> > mlt = {0.0, 0.0, 0.0};
        std::vector<std::complex<double> > mqt = {0.0, 0.0, 0.0};
        std::vector<std::complex<double> > mut = {0.0, 0.0, 0.0};

        std::vector<std::vector<std::complex<double> > > ad = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > adc = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > ae = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > aec = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > au = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > auc = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > Yd = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > Ydc = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > Ye = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > Yec = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > Yu = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        std::vector<std::vector<std::complex<double> > > Yuc = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    public:
        MSSM() = default;

        MSSM(std::unordered_map<std::string, std::complex<double> > params, double scale, bool loop);

        void updateParams(std::unordered_map<std::string, std::complex<double> > params);

        double getScale();

        void setScale(double scale);

        void loopContributions(bool loop);

        std::unordered_map<std::string, std::complex<double> > getParams();

        std::complex<double> cG();
        std::complex<double> cGt();
        std::complex<double> cH();
        std::complex<double> cHB();
        std::complex<double> cHBox();
        std::complex<double> cHBt();
        std::complex<double> cHD();
        std::complex<double> cHG();
        std::complex<double> cHGt();
        std::complex<double> cHW();
        std::complex<double> cHWB();
        std::complex<double> cHWt();
        std::complex<double> cHWtB();
        std::complex<double> cW();
        std::complex<double> cWt();
        std::complex<double> cdB(int i1, int i2);
        std::complex<double> cdG(int i1, int i2);
        std::complex<double> cdH(int i1, int i2);
        std::complex<double> cdW(int i1, int i2);
        std::complex<double> ceB(int i1, int i2);
        std::complex<double> ceH(int i1, int i2);
        std::complex<double> ceW(int i1, int i2);
        std::complex<double> cHd(int i1, int i2);
        std::complex<double> cHe(int i1, int i2);
        std::complex<double> cHl1(int i1, int i2);
        std::complex<double> cHl3(int i1, int i2);
        std::complex<double> cHq1(int i1, int i2);
        std::complex<double> cHq3(int i1, int i2);
        std::complex<double> cHu(int i1, int i2);
        std::complex<double> cHud(int i1, int i2);
        std::complex<double> cllHH(int i1, int i2);
        std::complex<double> cuB(int i1, int i2);
        std::complex<double> cuG(int i1, int i2);
        std::complex<double> cuH(int i1, int i2);
        std::complex<double> cuW(int i1, int i2);
        std::complex<double> cdd(int i1, int i2, int i3, int i4);
        std::complex<double> cduq(int i1, int i2, int i3, int i4);
        std::complex<double> cduu(int i1, int i2, int i3, int i4);
        std::complex<double> ced(int i1, int i2, int i3, int i4);
        std::complex<double> cee(int i1, int i2, int i3, int i4);
        std::complex<double> ceu(int i1, int i2, int i3, int i4);
        std::complex<double> cld(int i1, int i2, int i3, int i4);
        std::complex<double> cle(int i1, int i2, int i3, int i4);
        std::complex<double> cledq(int i1, int i2, int i3, int i4);
        std::complex<double> clequ1(int i1, int i2, int i3, int i4);
        std::complex<double> clequ3(int i1, int i2, int i3, int i4);
        std::complex<double> cll(int i1, int i2, int i3, int i4);
        std::complex<double> clq1(int i1, int i2, int i3, int i4);
        std::complex<double> clq3(int i1, int i2, int i3, int i4);
        std::complex<double> clu(int i1, int i2, int i3, int i4);
        std::complex<double> cqd1(int i1, int i2, int i3, int i4);
        std::complex<double> cqd8(int i1, int i2, int i3, int i4);
        std::complex<double> cqe(int i1, int i2, int i3, int i4);
        std::complex<double> cqq1(int i1, int i2, int i3, int i4);
        std::complex<double> cqq3(int i1, int i2, int i3, int i4);
        std::complex<double> cqqq(int i1, int i2, int i3, int i4);
        std::complex<double> cqqu(int i1, int i2, int i3, int i4);
        std::complex<double> cqu1(int i1, int i2, int i3, int i4);
        std::complex<double> cqu8(int i1, int i2, int i3, int i4);
        std::complex<double> cquqd1(int i1, int i2, int i3, int i4);
        std::complex<double> cquqd8(int i1, int i2, int i3, int i4);
        std::complex<double> cud1(int i1, int i2, int i3, int i4);
        std::complex<double> cud8(int i1, int i2, int i3, int i4);
        std::complex<double> cuu(int i1, int i2, int i3, int i4);

        std::map<std::string, std::complex<double> > batch_eval(const std::vector<Task>& tasks);
};

