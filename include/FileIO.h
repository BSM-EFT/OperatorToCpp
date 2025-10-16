/**
 * @file FileIO.h
 * @author Suraj Prakash
 * @date 2025-10-17
 * @brief A suite of utility functions to aid in reading input from and writing output to files
 */

#pragma once
#include <ios>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <utility>
#include <sstream>
#include <iomanip>
#include <fstream>

#define hb 0.006332574

/// enum to indicate the nature of the output with respect to the order of loop-expansion
enum class ORDER {
  TREE,
  FULL,
  SPLIT,
};

/// function to trim leading whitespaces in a string
std::string trim_left(std::string s);

/// function to trim trailing whitespaces in a string
std::string trim_right(std::string s);

/// function to trim leading as well as trailing whitespaces
std::string trim(std::string s);


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
 * @overload
 */
void read_params(std::string fname, std::unordered_map<std::string, double>& p_dict);

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
template<typename Model>
double eval_wc(Model m, std::string s, double mubarsq, double hbar) {
    double res{};

    std::pair<std::string,std::vector<int>> name_idx = split_name_idx(s);
    std::string name = std::get<0>(name_idx);
    std::vector<int> idx = std::get<1>(name_idx);

    switch (idx.size()) {
        case 0:
            res = m.fname_map_0f[name](mubarsq, hbar);
            break;
        case 2:
            res = m.fname_map_2f[name](idx[0]-1, idx[1]-1, mubarsq, hbar);
            break;
        case 4:
            res = m.fname_map_4f[name](idx[0]-1, idx[1]-1, idx[2]-1, idx[3]-1, mubarsq, hbar);
            break;
    }

    return res;
}

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
template<typename Model>
std::string create_row(Model& m, std::vector<std::string>& keys, std::vector<double>& vals, std::vector<std::string>& wc_names, double mubarsq, ORDER ord) {
    std::unordered_map<std::string, double> param_dict;
    for (int i = 0; i < keys.size(); ++i) param_dict.emplace(keys[i], vals[i]);
    m.updateParams(param_dict);

    std::ostringstream rowstream;
    rowstream << std::fixed << std::setprecision(2);
    for (double val: vals) rowstream << val << ",";
    rowstream << std::scientific << std::setprecision(1);

    if (ord == ORDER::TREE) {
        for (std::string wc: wc_names) rowstream << eval_wc(m, wc, mubarsq, 0.0) << ",";
    } else if (ord == ORDER::FULL) {
        for (std::string wc: wc_names) rowstream << eval_wc(m, wc, mubarsq, 0.006332574) << ",";
    } else {
        for (std::string wc: wc_names) rowstream << "[" << eval_wc(m, wc, mubarsq, 0.0) << "," << eval_wc(m, wc, mubarsq, 0.006332574) - eval_wc(m, wc, mubarsq, 0.0) << "],";
    }

    std::string row = rowstream.str();
    return row.substr(0, row.length()-1);
}

/**
 * Helper function to check if a file already exists. Prints a message to the user stating that the file will be overwritten.
 *
 * @param f Name of the file.
 */
void check_file_exists(std::string f);

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
template<typename Model>
void write_to_csv(std::string fname, Model& m, std::map<std::string, std::vector<double> > p_range_dict, std::vector<std::string> wc_names, double mubarsq, ORDER ord) {
    std::vector<std::string> keys;
    for(auto it = p_range_dict.begin(); it != p_range_dict.end(); ++it) keys.emplace_back(it->first);
    std::vector<std::vector<double> > p_combs = create_param_combs(p_range_dict);

    check_file_exists(fname);
    std::ofstream f1;
    f1.open(fname, std::ios::out | std::ios::app);
    std::ostringstream h_stream;
    for (std::string key: keys) h_stream << key << ",";
    for (std::string wc: wc_names) h_stream << wc << ",";

    std::string header = h_stream.str();
    f1 << header.substr(0, header.length()-1) << "\n";
    for (auto p_comb: p_combs) f1 << create_row(m, keys, p_comb, wc_names, mubarsq, ord) << "\n";
    f1.close();
}

/**
 * @overload
 * For the case when combinations of parameter values are already given instead of ranges for each parameter.
 */
template<typename Model>
void write_to_csv(std::string fname, Model& m, std::vector<std::vector<double> > p_combs, std::vector<std::string> keys, std::vector<std::string> wc_names, double mubarsq, ORDER ord) {
    check_file_exists(fname);
    std::ofstream f1;
    f1.open(fname, std::ios::out | std::ios::app);
    std::ostringstream h_stream;
    for (std::string key: keys) h_stream << key << ",";
    for (std::string wc: wc_names) h_stream << wc << ",";

    std::string header = h_stream.str();
    f1 << header.substr(0, header.length()-1) << "\n";
    for (auto p_comb: p_combs) f1 << create_row(m, keys, p_comb, wc_names, mubarsq, ord) << "\n";
    f1.close();
}

/**
 * Create a (.yaml) file that stores a set of combinations of specified parameters and evaluated coefficients for a single benchmark point.
 *
 * @param fname Name of the (.csv) file that will store the results.
 * @param m The physics model that defines the parameters as member variables and the coefficients as methods.
 * @param p_dict A map that stores (name, value) pairs for all model parameters.
 * @param keys A vector of names of the model parameters that are to be included in the output file.
 * @param wc_names A vector of coefficient names to be evaluated and added to the file.
 * @param mubarsq Square of the Mass/energy scale at which the evaluation occurs.
 * @param ord The order of the loop expansion, either tree, full or split (tree, loop) output is created.
 */
template<typename Model>
void write_to_yaml(std::string fname, Model& m, std::unordered_map<std::string, double> p_dict, std::vector<std::string> keys, std::vector<std::string> wc_names, double mubarsq, ORDER ord) {
    check_file_exists(fname);
    std::ofstream f1;
    f1.open(fname, std::ios::out | std::ios::app);

    f1 << std::fixed << std::setprecision(2);
    for (std::string key: keys) {
        f1 << key << ": " << p_dict[key] << "\n";
    }

    f1 << std::scientific << std::setprecision(1);
    if (ord == ORDER::TREE) {
        for (std::string wc: wc_names) f1 << wc << ": " << eval_wc(m, wc, mubarsq, 0.0) << "\n";
    } else if (ord == ORDER::FULL) {
        for (std::string wc: wc_names) f1 << wc << ": "  << eval_wc(m, wc, mubarsq, 0.006332574) << "\n";
    } else {
        for (std::string wc: wc_names) f1 << wc << ": "  << "[" << eval_wc(m, wc, mubarsq, 0.0) << "," << eval_wc(m, wc, mubarsq, 0.006332574) - eval_wc(m, wc, mubarsq, 0.0) << "]\n";
    }

    f1.close();
}
