#pragma once
#include <vector>
#include <string>
#include <unordered_map>

class MSSM {
    private:
        double cgamma = 0.0;
        double g1 = 0.0;
        double g2 = 0.0;
        double g3 = 0.0;
        double m1 = 0.0;
        double m2 = 0.0;
        double m3 = 0.0;
        double mPhi = 0.0;
        double muTilde = 0.0;

        std::vector<double> mdt = {0.0, 0.0, 0.0};
        std::vector<double> met = {0.0, 0.0, 0.0};
        std::vector<double> mlt = {0.0, 0.0, 0.0};
        std::vector<double> mqt = {0.0, 0.0, 0.0};
        std::vector<double> mut = {0.0, 0.0, 0.0};

        std::vector<std::vector<double> > ad = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
        std::vector<std::vector<double> > ae = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
        std::vector<std::vector<double> > au = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
        std::vector<std::vector<double> > yd = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
        std::vector<std::vector<double> > ye = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
        std::vector<std::vector<double> > yu = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};

    public:
        MSSM() = default;

        MSSM(std::unordered_map<std::string, double> params);

        void updateParams(std::unordered_map<std::string, double> params);

        void printParamNames();

        void printParams();

        double cllHH(int i1, int i2, double mubarsq, double hbar);
        double cG(double mubarsq, double hbar);
        double cW(double mubarsq, double hbar);
        double cGt(double mubarsq, double hbar);
        double cWt(double mubarsq, double hbar);
        double cH(double mubarsq, double hbar);
        double cHBox(double mubarsq, double hbar);
        double cHD(double mubarsq, double hbar);
        double cHG(double mubarsq, double hbar);
        double cHW(double mubarsq, double hbar);
        double cHB(double mubarsq, double hbar);
        double cHWB(double mubarsq, double hbar);
        double cHGt(double mubarsq, double hbar);
        double cHWt(double mubarsq, double hbar);
        double cHBt(double mubarsq, double hbar);
        double cHWtB(double mubarsq, double hbar);
        double ceH(int i1, int i2, double mubarsq, double hbar);
        double cuH(int i1, int i2, double mubarsq, double hbar);
        double cdH(int i1, int i2, double mubarsq, double hbar);
        double ceW(int i1, int i2, double mubarsq, double hbar);
        double ceB(int i1, int i2, double mubarsq, double hbar);
        double cuG(int i1, int i2, double mubarsq, double hbar);
        double cuW(int i1, int i2, double mubarsq, double hbar);
        double cuB(int i1, int i2, double mubarsq, double hbar);
        double cdG(int i1, int i2, double mubarsq, double hbar);
        double cdW(int i1, int i2, double mubarsq, double hbar);
        double cdB(int i1, int i2, double mubarsq, double hbar);
        double cHl1(int i1, int i2, double mubarsq, double hbar);
        double cHl3(int i1, int i2, double mubarsq, double hbar);
        double cHe(int i1, int i2, double mubarsq, double hbar);
        double cHq1(int i1, int i2, double mubarsq, double hbar);
        double cHq3(int i1, int i2, double mubarsq, double hbar);
        double cHu(int i1, int i2, double mubarsq, double hbar);
        double cHd(int i1, int i2, double mubarsq, double hbar);
        double cHud(int i1, int i2, double mubarsq, double hbar);
        double cll(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqq1(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqq3(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double clq1(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double clq3(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cee(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cuu(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cdd(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double ceu(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double ced(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cud1(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cud8(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cle(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double clu(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cld(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqe(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqu1(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqu8(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqd1(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqd8(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cledq(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cquqd1(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cquqd8(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double clequ1(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double clequ3(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cduq(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqqu(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cqqq(int i1, int i2, int i3, int i4, double mubarsq, double hbar);
        double cduu(int i1, int i2, int i3, int i4, double mubarsq, double hbar);

};

