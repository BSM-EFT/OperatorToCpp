/**
 * @file OperatorImport.h
 * @author Suraj Prakash
 * @date 2025-03-18
 * @brief header file correspondig to OperatorImport.cpp
 */

#include <vector>
#include <variant>
#include <tuple>
#include <functional>

typedef std::tuple<std::function<double(int, int, double)>, double> YF_tuple;

class LoopFunc {
    public:
        std::vector<std::variant<std::vector<double>, double> > masses;
        int code;
        double mubarsq;

        LoopFunc() = delete;
        LoopFunc(std::vector<std::variant<std::vector<double>, double> > list_of_masses, int code, double mubarsq);
};

double LF(std::vector<double> masses, int code, double mubarsq);

class MassPow {
    public:
        std::variant<std::vector<double>, double> mass;
        int exp;

        MassPow() = delete;
        MassPow(std::variant<std::vector<double>, double> mass, int exp);

};

double Eval(LoopFunc loopf, const std::vector<int>& idx);

double Eval(MassPow masspw, const std::vector<int>& idx);

double Eval(const std::vector<std::vector<double> >& matrix, const std::vector<int>& idx);

double Eval(YF_tuple x, std::vector<int> idx);

int dim(std::vector<double>& m);

int dim(double& m);

double apply(std::vector<double>& m, int i);

double apply(double& m, int i);

double exponentiate(std::vector<double> m, const std::vector<int> idx, int pw);

double exponentiate(double m, const std::vector<int> idx, int pw);

int maxRepIdx(const std::vector<std::vector<int> >& v1);

std::vector<std::vector<int> > idx_seqs(const std::vector<std::vector<int> >& index_seqs, const std::vector<int>& free_indices, const std::vector<int>& rep_idx_vals);

std::vector<std::vector<int> > idx_seqs(const std::vector<std::vector<int> >& index_seqs, const std::vector<int>& free_indices);

std::vector<std::vector<int> > cartesianProduct(int num_flavours, int num_idx);

double EinsSum(std::vector<std::variant<LoopFunc, MassPow, std::vector<std::vector<double> >, YF_tuple> > tensor_objs, std::vector<std::vector<int> > index_order, std::vector<int> free_indices);

int KronDelta(int a, int b);
