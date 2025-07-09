/**
 * @file FileIO.h
 * @author Suraj Prakash
 * @date 2025-07-09
 * @brief A suite of utility functions to aid in reading input from and writing output to files
 */

#include "MSSM.h"

using Model = MSSM;

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <utility>

/// enum to indicate the nature of the output with respect to the order of loop-expansion
enum class ORDER {
  TREE,
  FULL,
  SPLIT,
};

/// create a range from start to stop with step size delta.
std::vector<double> create_range(double start, double end, double delta);

/**
 * Read parameter values from a file and store them in maps.
 *
 * @param fname File whose contents are read (.yaml format is assumed).
 * @param p_dict An unordered map to store (name, value) pairs for parameters with fixed values.
 * @param p_range_dict A map to store (name, range) pairs for parameters whose min and max values are specified.
 */
void read_params(std::string fname, std::unordered_map<std::string, double>& p_dict, std::map<std::string, std::vector<double> >& p_range_dict);

/**
 * Create combinations of parameter values as cartesian products of parameter ranges.
 *
 * @param p_range_dict A map that stores (name, range) pairs for a set of parameters.
 * @return A vector of vectors, where each each inner vector contains a single combination of parameter values
 * following the order in which they appear as keys within p_range_dict.
 */
std::vector<std::vector<double> > create_param_combs(std::map<std::string, std::vector<double> >& p_range_dict);

/**
 * Create a vector of strings based on the "coefficients" listed in a .txt file.
 *
 * @param fname A file that stores coefficient names line by line..
 * @return A vector of strings corresponding to the names of individual coefficients.
 */
std::vector<std::string> read_wc_names(std::string fname);

/**
 * Separate the name and index information from the full coefficient name.
 *
 * @param full_name Coefficient name in the form "cABC" or "cABC_ijkl".
 * @return A tuple containing (name, index) pairs.
 */
std::pair<std::string,std::vector<int>> split_name_idx(std::string full_name);

/**
 * Evaluate the coefficient given its name and the high energy physics model.
 *
 * @param m The physics model that defines the coefficient as one of its methods.
 * @param s String containing the full name of the coefficient.
 * @param mubarsq square of the Mass/energy scale at which the evaluation occurs.
 * @return Numerical value of the coefficient.
 */
double eval_wc(Model m, std::string s, double mubarsq);

/**
 * Create a line containing values of independent parameters and evaluated coefficients.
 *
 * @param m The physics model that defines the parameters as member variables and the coefficients as methods.
 * @param keys A vector of parameter names.
 * @param vals A particular combination of parameter values corresponding to the keys.
 * @param wc_names A vector of coefficient names to be evaluated and added to the line.
 * @param mubarsq square of the Mass/energy scale at which the evaluation occurs.
 * @return String constituting a single line/row of parameter(s), coefficient(s) information.
 */
std::string create_row(Model& m, std::vector<std::string>& keys, std::vector<double>& vals, std::vector<std::string>& wc_names, double mubarsq);

/**
 * Create a (.csv) file that stores a set of combinations of independent parameters and evaluated coefficients for each combination.
 *
 * @param fname Name of the (.csv) file that will store the results.
 * @param m The physics model that defines the parameters as member variables and the coefficients as methods.
 * @param p_range_dict A map that stores (name, range) pairs for a set of parameters.
 * @param wc_names A vector of coefficient names to be evaluated and added to the file.
 * @param mubarsq Square of the Mass/energy scale at which the evaluation occurs.
 * @param ord The order of the loop expansion, either tree, full or split (tree, loop) output is created.
 */
void write_to_csv(std::string fname, Model& m, std::map<std::string, std::vector<double> > p_range_dict, std::vector<std::string> wc_names, double mubarsq, ORDER ord);
