/**
 * @file OperatorImport.cpp
 * @author Suraj Prakash
 * @date 2025-10-13
 * @brief Auxiliary classes and funtions to aid in the evaluation of expressions within the Wilson coefficient functions
 */

#include "OperatorImport.h"
#include "LF.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <variant>
#include <cassert>
#include <stdexcept>
#include <tuple>
#include <functional>
#include <string>

using std::vector;
using std::transform;
using std::accumulate;
using std::variant;
using std::invalid_argument;
using std::string;

LoopFunc::LoopFunc(vector<variant<vector<double>, double> > list_of_masses, int code, double mubarsq) {
  for (variant<vector<double>, double> mass_vec: list_of_masses) {
    this->masses.emplace_back(mass_vec);
  }
  this->code = code;
  this->mubarsq = mubarsq;
}

MassPow::MassPow(variant<vector<double>, double> mass, int exp) {
    this->mass = mass;
    this->exp = exp;
}

// Eval function that operates on a matrix and a vector containing indices {i, j}
double Eval(const vector<vector<double> >& matrix, const vector<int>& idx) {
    assert(idx.size()==2);
    return matrix[idx[0]][idx[1]];
}

// overloaded Eval function for the special case of Yukawa functions
typedef std::tuple<std::function<double(int, int, double)>, double> YF_tuple;

double Eval(YF_tuple x, vector<int> idx) {
   std::function<double(int, int, double)> f = std::get<0>(x);
   double mubarsq = std::get<1>(x);
   return f(idx[0], idx[1], mubarsq);
}

// overloaded Eval function that returns the appropriate loop function for given combination of masses and exponents
double Eval(LoopFunc loopf, const vector<int>& idx) {
    vector<double> mass_arg;
    int i = 0;

    for (auto mass : loopf.masses) {
        std::visit([&mass_arg, &idx, &i](auto m){
            if (dim(m)==0) { mass_arg.emplace_back(apply(m,0)); }
            else {
                mass_arg.emplace_back(apply(m,idx[i]));
                i++;
            }
        }, mass);
    }
    return LF(mass_arg, loopf.code, loopf.mubarsq);
}

// overloaded Eval function that operates on a MassPow object and a vector containing zero or one element
double Eval(MassPow masspw, const vector<int>& idx) {
    variant<vector<double>, double> mass = masspw.mass;
    int pw = masspw.exp;
    double res = 1;

    std::visit([&pw, &idx, &res](auto m){ res = exponentiate(m, idx, pw); }, mass);
    return res;
}

// helper to helper functions
int dim(vector<double>& m) { return 1; }

int dim(double& m) { return 0; }

double apply(vector<double>& m, int i) { return m[i]; }

double apply(double& m, int i) {
    assert(i==0);
    return m;
}

double exponentiate(vector<double> m, const vector<int> idx, int pw) {
    assert(idx.size()==1);
    return pow(m[idx[0]], pw);
}

double exponentiate(double m, const vector<int> idx, int pw) {
    assert(idx.size()==0);
    return pow(m, pw);
}

// function to identify the number of repeated indices
int maxRepIdx(const vector<vector<int> >& v1) {
    vector<int> v2;

    for(vector<int> v : v1){
        if (v.size()==0) continue;
        else {
            transform(v.begin(), v.end(), v.begin(), [](int i){return i>10 ? 0 : i%10;});
            v2.emplace_back(*max_element(v.begin(), v.end()));
        }
    };

    if (v2.size()==0) return 0;
    else return *max_element(v2.begin(), v2.end());
}

// building a vector of vector of index values while keeping track of the free-indices appropriately.
vector<vector<int> > idx_seqs(const vector<vector<int> >& index_seqs, const vector<int>& free_indices, const vector<int>& rep_idx_vals){
    vector<vector<int> > addrs;
    int free_ind_pos = 0;

    for(vector<int> v : index_seqs) {
        vector<int> dest;
        for(int i : v) {
            if (i>10 && (free_indices.size() > free_ind_pos)) {
                dest.emplace_back(free_indices[free_ind_pos++]);
            } else if (i > 10) {
                throw invalid_argument {"Mismatch in the number of free indices across arguemnts of idx_seqs."};
            } else dest.emplace_back(rep_idx_vals[i-1]);
        }
        addrs.emplace_back(dest);
    }
    return addrs;
}

// overloaded version of idx_seqs for the case where there are no repeated indices
vector<vector<int> > idx_seqs(const vector<vector<int> >& index_seqs, const vector<int>& free_indices){
    vector<vector<int> > addrs;
    int free_ind_pos = 0;

    for(vector<int> v : index_seqs) {
        vector<int> dest;
        for(int i : v) {
            if (i>10 && (free_indices.size() > free_ind_pos)) dest.emplace_back(free_indices[free_ind_pos++]);
            else if (i > 10) throw invalid_argument {"Mismatch in the number of free indices across arguemnts of idx_seqs."};
        }
        addrs.emplace_back(dest);
    }
    return addrs;
}

// generating all possible combinations of values assumed by the repeated indices
vector<vector<int> > cartesianProduct(int num_flavours, int num_idx) {
    assert(num_idx > 0);
    int N = pow(num_flavours, num_idx);
    vector<vector<int> > allProducts;

    for(int i = 0; i < N; i++) {
        vector<int> vec(num_idx);
        int num = i;
        for(int j = num_idx-1; j >= 0; j--) {
            vec[j] = num % num_flavours;
            num = num / num_flavours;
        }
        allProducts.emplace_back(vec);
    }
    return allProducts;
}

// Einstein summation function for a variant containing (a loop function, a mass raised to some power, matrices or the Yukawa functions)
// and a specified ordering of repeated & free indices
double EinsSum(vector<variant<LoopFunc, MassPow, vector<vector<double> >, YF_tuple> > tensor_objs, vector<vector<int> > index_order, vector<int> free_indices) {
    int num_flavours = 3;   // hardcoded for now, should never be 0 or negative
    int num_idx = maxRepIdx(index_order);

    double sum{};

    if (num_idx == 0) { // if there are no repeated indices, we only need to evaluate once, no looping necessary
        vector<vector<int> > segregated_seqs = idx_seqs(index_order, free_indices);
        vector<double> res1;

        for (int k = 0; k < tensor_objs.size(); k++) {
            std::visit([&res1, &segregated_seqs, &k](auto obj){
                res1.emplace_back(Eval(obj, segregated_seqs[k]));
                }, tensor_objs[k]);
        }
        sum = accumulate(res1.begin(), res1.end(), 1.0, std::multiplies<double>());

    } else {
        vector<vector<int> > cprod = cartesianProduct(num_flavours, num_idx);

        for (vector<int> seq : cprod) {
            vector<vector<int> > segregated_seqs = idx_seqs(index_order, free_indices, seq);
            vector<double> res1;
            for (int k = 0; k < tensor_objs.size(); k++)
                std::visit([&res1, &segregated_seqs, &k](auto obj){res1.emplace_back(Eval(obj, segregated_seqs[k]));}, tensor_objs[k]);

            double res2 = accumulate(res1.begin(), res1.end(), 1.0, std::multiplies<double>());
            sum += res2;
        }
    }
    return sum;
}

// Kronecker delta function
int KronDelta(int a, int b) {
  if (a==b) return 1;
  else return 0;
}
