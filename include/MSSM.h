
#include <vector>
#include <string>
#include <map>
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

        MSSM(std::map<std::string, double> params);

        void updateParams(std::map<std::string, double> params);

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


};
