#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>


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

    std::unordered_map<std::string, std::function<double(double,double)>> fname_map_0f = {
        {"cG", [this](double msq, double hbar){ return this ->cG(msq, hbar); } },
        {"cW", [this](double msq, double hbar){ return this ->cW(msq, hbar); } },
        {"cGt", [this](double msq, double hbar){ return this ->cGt(msq, hbar); } },
        {"cWt", [this](double msq, double hbar){ return this ->cWt(msq, hbar); } },
        {"cH", [this](double msq, double hbar){ return this ->cH(msq, hbar); } },
        {"cHBox", [this](double msq, double hbar){ return this ->cHBox(msq, hbar); } },
        {"cHD", [this](double msq, double hbar){ return this ->cHD(msq, hbar); } },
        {"cHG", [this](double msq, double hbar){ return this ->cHG(msq, hbar); } },
        {"cHW", [this](double msq, double hbar){ return this ->cHW(msq, hbar); } },
        {"cHB", [this](double msq, double hbar){ return this ->cHB(msq, hbar); } },
        {"cHWB", [this](double msq, double hbar){ return this ->cHWB(msq, hbar); } },
        {"cHGt", [this](double msq, double hbar){ return this ->cHGt(msq, hbar); } },
        {"cHWt", [this](double msq, double hbar){ return this ->cHWt(msq, hbar); } },
        {"cHBt", [this](double msq, double hbar){ return this ->cHBt(msq, hbar); } },
        {"cHWtB", [this](double msq, double hbar){ return this ->cHWtB(msq, hbar); } }
    };

    std::unordered_map<std::string, std::function<double(int,int,double,double)>> fname_map_2f = {
        {"cllHH", [this](int k, int l, double msq, double hbar){ return this ->cllHH(k, l, msq, hbar); } },
        {"ceH", [this](int k, int l, double msq, double hbar){ return this ->ceH(k, l, msq, hbar); } },
        {"cuH", [this](int k, int l, double msq, double hbar){ return this ->cuH(k, l, msq, hbar); } },
        {"cdH", [this](int k, int l, double msq, double hbar){ return this ->cdH(k, l, msq, hbar); } },
        {"ceW", [this](int k, int l, double msq, double hbar){ return this ->ceW(k, l, msq, hbar); } },
        {"ceB", [this](int k, int l, double msq, double hbar){ return this ->ceB(k, l, msq, hbar); } },
        {"cuG", [this](int k, int l, double msq, double hbar){ return this ->cuG(k, l, msq, hbar); } },
        {"cuW", [this](int k, int l, double msq, double hbar){ return this ->cuW(k, l, msq, hbar); } },
        {"cuB", [this](int k, int l, double msq, double hbar){ return this ->cuB(k, l, msq, hbar); } },
        {"cdG", [this](int k, int l, double msq, double hbar){ return this ->cdG(k, l, msq, hbar); } },
        {"cdW", [this](int k, int l, double msq, double hbar){ return this ->cdW(k, l, msq, hbar); } },
        {"cdB", [this](int k, int l, double msq, double hbar){ return this ->cdB(k, l, msq, hbar); } },
        {"cHl1", [this](int k, int l, double msq, double hbar){ return this ->cHl1(k, l, msq, hbar); } },
        {"cHl3", [this](int k, int l, double msq, double hbar){ return this ->cHl3(k, l, msq, hbar); } },
        {"cHe", [this](int k, int l, double msq, double hbar){ return this ->cHe(k, l, msq, hbar); } },
        {"cHq1", [this](int k, int l, double msq, double hbar){ return this ->cHq1(k, l, msq, hbar); } },
        {"cHq3", [this](int k, int l, double msq, double hbar){ return this ->cHq3(k, l, msq, hbar); } },
        {"cHu", [this](int k, int l, double msq, double hbar){ return this ->cHu(k, l, msq, hbar); } },
        {"cHd", [this](int k, int l, double msq, double hbar){ return this ->cHd(k, l, msq, hbar); } },
        {"cHud", [this](int k, int l, double msq, double hbar){ return this ->cHud(k, l, msq, hbar); } }
    };

    std::unordered_map<std::string, std::function<double(int,int,int,int,double,double)>> fname_map_4f = {
        {"cll", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cll(i, j, k, l, msq, hbar); } },
        {"cqq1", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqq1(i, j, k, l, msq, hbar); } },
        {"cqq3", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqq3(i, j, k, l, msq, hbar); } },
        {"clq1", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->clq1(i, j, k, l, msq, hbar); } },
        {"clq3", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->clq3(i, j, k, l, msq, hbar); } },
        {"cee", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cee(i, j, k, l, msq, hbar); } },
        {"cuu", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cuu(i, j, k, l, msq, hbar); } },
        {"cdd", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cdd(i, j, k, l, msq, hbar); } },
        {"ceu", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->ceu(i, j, k, l, msq, hbar); } },
        {"ced", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->ced(i, j, k, l, msq, hbar); } },
        {"cud1", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cud1(i, j, k, l, msq, hbar); } },
        {"cud8", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cud8(i, j, k, l, msq, hbar); } },
        {"cle", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cle(i, j, k, l, msq, hbar); } },
        {"clu", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->clu(i, j, k, l, msq, hbar); } },
        {"cld", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cld(i, j, k, l, msq, hbar); } },
        {"cqe", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqe(i, j, k, l, msq, hbar); } },
        {"cqu1", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqu1(i, j, k, l, msq, hbar); } },
        {"cqu8", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqu8(i, j, k, l, msq, hbar); } },
        {"cqd1", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqd1(i, j, k, l, msq, hbar); } },
        {"cqd8", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqd8(i, j, k, l, msq, hbar); } },
        {"cledq", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cledq(i, j, k, l, msq, hbar); } },
        {"cquqd1", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cquqd1(i, j, k, l, msq, hbar); } },
        {"cquqd8", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cquqd8(i, j, k, l, msq, hbar); } },
        {"clequ1", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->clequ1(i, j, k, l, msq, hbar); } },
        {"clequ3", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->clequ3(i, j, k, l, msq, hbar); } },
        {"cduq", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cduq(i, j, k, l, msq, hbar); } },
        {"cqqu", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqqu(i, j, k, l, msq, hbar); } },
        {"cqqq", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cqqq(i, j, k, l, msq, hbar); } },
        {"cduu", [this](int i, int j, int k, int l, double msq, double hbar){ return this ->cduu(i, j, k, l, msq, hbar); } }
    };

};

