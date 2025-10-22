/**
 * @file OperatorImport.h
 * @author Suraj Prakash
 * @date 2025-10-17
 * @brief Auxiliary classes and funtions to aid in the evaluation of expressions within the Wilson coefficient functions
 */

#pragma once
#include <vector>
#include <string>
#include <variant>
#include <tuple>
#include <utility>
#include <functional>

/// Type alias declaration for a (function object, double) tuple.
typedef std::tuple<std::function<double(int, int, double)>, double> YF_tuple;

/// Implemntation of a loop function class to faciliate calling and evaluating LoopFunc functors.
class LoopFunc {
    public:
        std::vector<std::variant<std::vector<double>, double> > masses;
        int code;
        double mubarsq;

        LoopFunc() = delete;
        LoopFunc(std::vector<std::variant<std::vector<double>, double> > list_of_masses, int code, double mubarsq);
};

/// Implemntation of a class to faciliate calling and evaluating MassPow functors.
class MassPow {
    public:
        std::variant<std::vector<double>, double> mass;
        int exp;

        MassPow() = delete;
        MassPow(std::variant<std::vector<double>, double> mass, int exp);

};

/**
 * Evaluate a LoopFunc object given a set of indices.
 *
 * @param loopf A LoopFunc object with a specific set of masses, function variety and mass scale for evaluation.
 * @param idx A vector of index values for the 1d "vector" masses within the LoopFunc mass list.
 * @return The evaluated numerical value.
 */
double Eval(LoopFunc loopf, const std::vector<int>& idx);

/**
 * Evaluate a MassPow object given a set of indices.
 *
 * @overload
 */
double Eval(MassPow masspw, const std::vector<int>& idx);

/**
 * Evaluate a 2d vector (matrix) given a set of indices.
 *
 * @overload
 */
double Eval(const std::vector<std::vector<double> >& matrix, const std::vector<int>& idx);

/**
 * Evaluate a (function object, mass scale) tuple given a set of indices.
 *
 * @overload
 */
double Eval(YF_tuple x, std::vector<int> idx);

/// return the dimension of a 1d vector -> necessary for std::visit calls.
int dim(std::vector<double>& m);

/// return the dimension of a 0d vector -> necessary for std::visit calls.
int dim(double& m);

/// subscript a 1d vector with a specified index -> necessary for std::visit calls.
double apply(std::vector<double>& m, int i);

/// overloaded apply function for double valued masses, simply returns the original value.
double apply(double& m, int i);

/// exponentiate a specified element of a 1d vector by a specified power.
double exponentiate(std::vector<double> m, const std::vector<int> idx, int pw);

/// overloaded exponentiate function for double valued masses, the index information is discarded.
double exponentiate(double m, const std::vector<int> idx, int pw);

/// identify the maximum number of unique indices given a vector of vector of indices.
int maxRepIdx(const std::vector<std::vector<int> >& v1);

/**
 * Create a sequence of index values based on a given permutation of indices, along with values for free and repeated indices.
 *
 * @param index_seqs The given permutation of free and repeated indices stored as a vector of vectors.
 * @param free_indices A vector storing the values of each free index.
 * @param rep_idx_vals A vector containing the current value of each repeated index.
 * @return A vector of vectors containing values of each index, repeated or free.
*/
std::vector<std::vector<int> > idx_seqs(const std::vector<std::vector<int> >& index_seqs, const std::vector<int>& free_indices, const std::vector<int>& rep_idx_vals);

/**
 * Create a sequence of index values based on a given permutation of indices, when there are no repeated indices.
 *
 * @overload
*/
std::vector<std::vector<int> > idx_seqs(const std::vector<std::vector<int> >& index_seqs, const std::vector<int>& free_indices);

/// Create a vector of all possible cartesian products of index values given the number of flavours and number of indices.
std::vector<std::vector<int> > cartesianProduct(int num_flavours, int num_idx);

/**
 * Evaluate the Einstein sum over a product of indical objects with a given permutation of repeated index contraction.
 *
 * @param tensor_objs A vector of variants, i.e. object (with indices) of LoopFunc, MassPow, matrix or YF_tuple types.
 * @param index_order The given permutation of free and repeated indices stored as a vector of vectors.
 * @param free_indices A vector storing the values of each free index.
 * @return The evaluated sum.
*/
double EinsSum(std::vector<std::variant<LoopFunc, MassPow, std::vector<std::vector<double> >, YF_tuple> > tensor_objs, std::vector<std::vector<int> > index_order, std::vector<int> free_indices);

/// Evaluate the kronecker delta function for specified indices, returns 1 for a == b, 0 otherwise.
int KronDelta(int a, int b);
