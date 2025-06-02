
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>
#define hbar 0.006332574


class MSSM {
    private:
        double cgamma = 0.0;
        double g1 = 0.0;
        double g2 = 0.0;
        double g3 = 0.0;
        double m1 = 0.0;
        double m2 = 0.0;
        double m3 = 0.0;
        double mHsq = 0.0;
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

        double cllHH(int i1, int i2, double mubarsq);
        double cG(double mubarsq);
        double cW(double mubarsq);
        double cGt(double mubarsq);
        double cWt(double mubarsq);
        double cH(double mubarsq);
        double cHBox(double mubarsq);
        double cHD(double mubarsq);
        double cHG(double mubarsq);
        double cHW(double mubarsq);
        double cHB(double mubarsq);
        double cHWB(double mubarsq);
        double cHGt(double mubarsq);
        double cHWt(double mubarsq);
        double cHBt(double mubarsq);
        double cHWtB(double mubarsq);
        double ceH(int i1, int i2, double mubarsq);
        double cuH(int i1, int i2, double mubarsq);
        double cdH(int i1, int i2, double mubarsq);
        double ceW(int i1, int i2, double mubarsq);
        double ceB(int i1, int i2, double mubarsq);
        double cuG(int i1, int i2, double mubarsq);
        double cuW(int i1, int i2, double mubarsq);
        double cuB(int i1, int i2, double mubarsq);
        double cdG(int i1, int i2, double mubarsq);
        double cdW(int i1, int i2, double mubarsq);
        double cdB(int i1, int i2, double mubarsq);
        double cHl1(int i1, int i2, double mubarsq);
        double cHl3(int i1, int i2, double mubarsq);
        double cHe(int i1, int i2, double mubarsq);
        double cHq1(int i1, int i2, double mubarsq);
        double cHq3(int i1, int i2, double mubarsq);
        double cHu(int i1, int i2, double mubarsq);
        double cHd(int i1, int i2, double mubarsq);
        double cHud(int i1, int i2, double mubarsq);
        double cll(int i1, int i2, int i3, int i4, double mubarsq);
        double cqq1(int i1, int i2, int i3, int i4, double mubarsq);
        double cqq3(int i1, int i2, int i3, int i4, double mubarsq);
        double clq1(int i1, int i2, int i3, int i4, double mubarsq);
        double clq3(int i1, int i2, int i3, int i4, double mubarsq);
        double cee(int i1, int i2, int i3, int i4, double mubarsq);
        double cuu(int i1, int i2, int i3, int i4, double mubarsq);
        double cdd(int i1, int i2, int i3, int i4, double mubarsq);
        double ceu(int i1, int i2, int i3, int i4, double mubarsq);
        double ced(int i1, int i2, int i3, int i4, double mubarsq);
        double cud1(int i1, int i2, int i3, int i4, double mubarsq);
        double cud8(int i1, int i2, int i3, int i4, double mubarsq);
        double cle(int i1, int i2, int i3, int i4, double mubarsq);
        double clu(int i1, int i2, int i3, int i4, double mubarsq);
        double cld(int i1, int i2, int i3, int i4, double mubarsq);
        double cqe(int i1, int i2, int i3, int i4, double mubarsq);
        double cqu1(int i1, int i2, int i3, int i4, double mubarsq);
        double cqu8(int i1, int i2, int i3, int i4, double mubarsq);
        double cqd1(int i1, int i2, int i3, int i4, double mubarsq);
        double cqd8(int i1, int i2, int i3, int i4, double mubarsq);
        double cledq(int i1, int i2, int i3, int i4, double mubarsq);
        double cquqd1(int i1, int i2, int i3, int i4, double mubarsq);
        double cquqd8(int i1, int i2, int i3, int i4, double mubarsq);
        double clequ1(int i1, int i2, int i3, int i4, double mubarsq);
        double clequ3(int i1, int i2, int i3, int i4, double mubarsq);
        double cduq(int i1, int i2, int i3, int i4, double mubarsq);
        double cqqu(int i1, int i2, int i3, int i4, double mubarsq);
        double cqqq(int i1, int i2, int i3, int i4, double mubarsq);
        double cduu(int i1, int i2, int i3, int i4, double mubarsq);

    std::unordered_map<std::string, std::function<double(double)>> fname_map_0f = {
        {"cG", [this](double msq){ return this ->cG(msq); } },
        {"cW", [this](double msq){ return this ->cW(msq); } },
        {"cGt", [this](double msq){ return this ->cGt(msq); } },
        {"cWt", [this](double msq){ return this ->cWt(msq); } },
        {"cH", [this](double msq){ return this ->cH(msq); } },
        {"cHBox", [this](double msq){ return this ->cHBox(msq); } },
        {"cHD", [this](double msq){ return this ->cHD(msq); } },
        {"cHG", [this](double msq){ return this ->cHG(msq); } },
        {"cHW", [this](double msq){ return this ->cHW(msq); } },
        {"cHB", [this](double msq){ return this ->cHB(msq); } },
        {"cHWB", [this](double msq){ return this ->cHWB(msq); } },
        {"cHGt", [this](double msq){ return this ->cHGt(msq); } },
        {"cHWt", [this](double msq){ return this ->cHWt(msq); } },
        {"cHBt", [this](double msq){ return this ->cHBt(msq); } },
        {"cHWtB", [this](double msq){ return this ->cHWtB(msq); } }
    };

    std::unordered_map<std::string, std::function<double(int,int,double)>> fname_map_2f = {
        {"cllHH", [this](int k, int l, double msq){ return this ->cllHH(k, l, msq); } },
        {"ceH", [this](int k, int l, double msq){ return this ->ceH(k, l, msq); } },
        {"cuH", [this](int k, int l, double msq){ return this ->cuH(k, l, msq); } },
        {"cdH", [this](int k, int l, double msq){ return this ->cdH(k, l, msq); } },
        {"ceW", [this](int k, int l, double msq){ return this ->ceW(k, l, msq); } },
        {"ceB", [this](int k, int l, double msq){ return this ->ceB(k, l, msq); } },
        {"cuG", [this](int k, int l, double msq){ return this ->cuG(k, l, msq); } },
        {"cuW", [this](int k, int l, double msq){ return this ->cuW(k, l, msq); } },
        {"cuB", [this](int k, int l, double msq){ return this ->cuB(k, l, msq); } },
        {"cdG", [this](int k, int l, double msq){ return this ->cdG(k, l, msq); } },
        {"cdW", [this](int k, int l, double msq){ return this ->cdW(k, l, msq); } },
        {"cdB", [this](int k, int l, double msq){ return this ->cdB(k, l, msq); } },
        {"cHl1", [this](int k, int l, double msq){ return this ->cHl1(k, l, msq); } },
        {"cHl3", [this](int k, int l, double msq){ return this ->cHl3(k, l, msq); } },
        {"cHe", [this](int k, int l, double msq){ return this ->cHe(k, l, msq); } },
        {"cHq1", [this](int k, int l, double msq){ return this ->cHq1(k, l, msq); } },
        {"cHq3", [this](int k, int l, double msq){ return this ->cHq3(k, l, msq); } },
        {"cHu", [this](int k, int l, double msq){ return this ->cHu(k, l, msq); } },
        {"cHd", [this](int k, int l, double msq){ return this ->cHd(k, l, msq); } },
        {"cHud", [this](int k, int l, double msq){ return this ->cHud(k, l, msq); } }
    };

    std::unordered_map<std::string, std::function<double(int,int,int,int,double)>> fname_map_4f = {
        {"cll", [this](int i, int j, int k, int l, double msq){ return this ->cll(i, j, k, l, msq); } },
        {"cqq1", [this](int i, int j, int k, int l, double msq){ return this ->cqq1(i, j, k, l, msq); } },
        {"cqq3", [this](int i, int j, int k, int l, double msq){ return this ->cqq3(i, j, k, l, msq); } },
        {"clq1", [this](int i, int j, int k, int l, double msq){ return this ->clq1(i, j, k, l, msq); } },
        {"clq3", [this](int i, int j, int k, int l, double msq){ return this ->clq3(i, j, k, l, msq); } },
        {"cee", [this](int i, int j, int k, int l, double msq){ return this ->cee(i, j, k, l, msq); } },
        {"cuu", [this](int i, int j, int k, int l, double msq){ return this ->cuu(i, j, k, l, msq); } },
        {"cdd", [this](int i, int j, int k, int l, double msq){ return this ->cdd(i, j, k, l, msq); } },
        {"ceu", [this](int i, int j, int k, int l, double msq){ return this ->ceu(i, j, k, l, msq); } },
        {"ced", [this](int i, int j, int k, int l, double msq){ return this ->ced(i, j, k, l, msq); } },
        {"cud1", [this](int i, int j, int k, int l, double msq){ return this ->cud1(i, j, k, l, msq); } },
        {"cud8", [this](int i, int j, int k, int l, double msq){ return this ->cud8(i, j, k, l, msq); } },
        {"cle", [this](int i, int j, int k, int l, double msq){ return this ->cle(i, j, k, l, msq); } },
        {"clu", [this](int i, int j, int k, int l, double msq){ return this ->clu(i, j, k, l, msq); } },
        {"cld", [this](int i, int j, int k, int l, double msq){ return this ->cld(i, j, k, l, msq); } },
        {"cqe", [this](int i, int j, int k, int l, double msq){ return this ->cqe(i, j, k, l, msq); } },
        {"cqu1", [this](int i, int j, int k, int l, double msq){ return this ->cqu1(i, j, k, l, msq); } },
        {"cqu8", [this](int i, int j, int k, int l, double msq){ return this ->cqu8(i, j, k, l, msq); } },
        {"cqd1", [this](int i, int j, int k, int l, double msq){ return this ->cqd1(i, j, k, l, msq); } },
        {"cqd8", [this](int i, int j, int k, int l, double msq){ return this ->cqd8(i, j, k, l, msq); } },
        {"cledq", [this](int i, int j, int k, int l, double msq){ return this ->cledq(i, j, k, l, msq); } },
        {"cquqd1", [this](int i, int j, int k, int l, double msq){ return this ->cquqd1(i, j, k, l, msq); } },
        {"cquqd8", [this](int i, int j, int k, int l, double msq){ return this ->cquqd8(i, j, k, l, msq); } },
        {"clequ1", [this](int i, int j, int k, int l, double msq){ return this ->clequ1(i, j, k, l, msq); } },
        {"clequ3", [this](int i, int j, int k, int l, double msq){ return this ->clequ3(i, j, k, l, msq); } },
        {"cduq", [this](int i, int j, int k, int l, double msq){ return this ->cduq(i, j, k, l, msq); } },
        {"cqqu", [this](int i, int j, int k, int l, double msq){ return this ->cqqu(i, j, k, l, msq); } },
        {"cqqq", [this](int i, int j, int k, int l, double msq){ return this ->cqqq(i, j, k, l, msq); } },
        {"cduu", [this](int i, int j, int k, int l, double msq){ return this ->cduu(i, j, k, l, msq); } }
    };

};

