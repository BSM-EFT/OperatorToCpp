/**
 * @file write_to_files.cpp
 * @author Suraj Prakash
 * @date 2025-06-02
 * @brief Code for writing WC values to file using MSSM.h and MSSM.cpp
 */

#include "MSSM.h"
#include "OperatorImport.h"
#include <ios>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>

#define mbarsq 1.0
#define step 0.1

using std::vector;
using std::unordered_map;
using std::map;
using std::string;
using std::ostringstream;
using std::setprecision;
using std::ofstream;
using std::ifstream;
using std::ios;

using Model = MSSM;

vector<double> create_range(double start, double end, double delta) {
    vector<double> vec;
    for (int i = 0; i < (end-start)/delta; i++)
        vec.emplace_back(start + i*delta);

    return vec;
}

void read_params(string fname, unordered_map<string, double>& p_dict, map<string, vector<double> >& p_range_dict) {
    ifstream params_file(fname);
    string line, p_name, rest;
    double val, start, end;

    while (std::getline(params_file, line)) {
        int sep = line.find(":");
        p_name = line.substr(0,sep);
        rest = line.substr(sep+1,line.length() - p_name.length());

        int b_open = rest.find("[");
        if (b_open == string::npos) {
            val = std::stod(rest.substr(1, rest.length()));
            p_dict.emplace(p_name, val);
        } else {
            int b_close = rest.find("]");
            int comma = rest.find(",");
            start = std::stod(rest.substr(b_open + 1, comma - b_open - 1));
            end = std::stod(rest.substr(comma + 2, b_close - b_open - 2));
            vector<double> val_vec = create_range(start, end, step);
            p_range_dict.emplace(p_name, val_vec);
        }
    }
}

vector<vector<double> > create_param_combs(map<string, vector<double> >& p_range_dict) {
    vector<vector<double>> val_combs, param_combs;
    vector<vector<int> > idx_combs = cartesianProduct(p_range_dict.begin()->second.size(), p_range_dict.size());

    for(auto it = p_range_dict.begin(); it != p_range_dict.end(); ++it) val_combs.emplace_back(it->second);

    for (vector<int> comb : idx_combs) {
        vector<double> par_comb;
        for (int i = 0; i < comb.size(); i++) {
            double val = val_combs[i][comb[i]];
            par_comb.emplace_back(val);
        }
        param_combs.emplace_back(par_comb);
    }
    return param_combs;
}

vector<string> read_wc_names(string fname) {
    ifstream wcs_file(fname);
    vector<string> wc_names;

    string line;
    while (std::getline(wcs_file, line)) wc_names.emplace_back(line);

    return wc_names;
}

string create_row(Model& m, vector<string>& keys, vector<double>& vals, vector<string>& wc_names) {
    unordered_map<string, double> param_dict;
    for (int i = 0; i < keys.size(); ++i) param_dict.emplace(keys[i], vals[i]);
    m.updateParams(param_dict);

    ostringstream rowstream;
    rowstream << std::fixed << setprecision(2);
    for (double val: vals) rowstream << val << ",";
    rowstream << std::scientific << setprecision(5);
    for (string wc: wc_names) rowstream << eval_wc(m, wc, mbarsq) << ",";

    string row = rowstream.str();
    return row.substr(0, row.length()-1);
}

void write_to_csv(Model& m, map<string, vector<double> > p_range_dict, vector<string> wc_names) {
    vector<string> keys;
    for(auto it = p_range_dict.begin(); it != p_range_dict.end(); ++it) keys.emplace_back(it->first);
    vector<vector<double> > p_combs = create_param_combs(p_range_dict);

    string fname = "./plots/data.csv";
    ofstream f1;
    f1.open(fname, ios::out | ios::app);
    ostringstream h_stream;
    for (string key: keys) h_stream << key << ",";
    for (string wc: wc_names) h_stream << wc << ",";

    string header = h_stream.str();
    f1 << header.substr(0, header.length()-1) << "\n";
    for (auto p_comb: p_combs) f1 << create_row(m, keys, p_comb, wc_names) << "\n";
    f1.close();
}

int main() {

    // create an instance of the model
    Model sb_model;

    // dictionaries to store fixed and variable value parameters
    unordered_map<string, double> par_dict;
    map<string, vector<double> > par_range_dict;

    // read parameter-values and wc-names as input from files
    read_params("params.yaml", par_dict, par_range_dict);
    vector<string> wcs = read_wc_names("coeffs.txt");

    // update the parameters of the model based on the fixed-valued inputs
    sb_model.updateParams(par_dict);

    // generate results and store them in a data.csv file
    write_to_csv(sb_model, par_range_dict, wcs);

    return 0;
}
