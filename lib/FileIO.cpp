/**
 * @file FileIO.cpp
 * @author Suraj Prakash
 * @date 2025-07-09
 * @brief A suite of utility functions to aid in reading input from and writing output to files
 */

#include "FileIO.h"
#include "OperatorImport.h"

using Model = MSSM;

#include <ios>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <utility>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <filesystem>

#define step 0.1
#define hb 0.006332574

using std::vector;
using std::unordered_map;
using std::map;
using std::string;
using std::ostringstream;
using std::setprecision;
using std::ofstream;
using std::ifstream;
using std::ios;


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

std::pair<string,vector<int>> split_name_idx(string full_name) {
    string name, idxStr;
    vector<int> idx;

    int sep = full_name.find("_");
    if (sep == string::npos) {
        name = full_name;
    } else {
        name = full_name.substr(0,sep);
        idxStr = full_name.substr(sep+1,full_name.length() - name.length());
        std::reverse(idxStr.begin(), idxStr.end());
        int num = std::stoi(idxStr);
        for (int i = 0; i < idxStr.length(); i++) {
            idx.emplace_back(num%10);
            num /= 10;
        }
    }

    return std::make_pair(name, idx);
}


double eval_wc(Model m, string s, double mubarsq, double hbar) {
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

string create_row(Model& m, vector<string>& keys, vector<double>& vals, vector<string>& wc_names, double mubarsq, ORDER ord) {
    unordered_map<string, double> param_dict;
    for (int i = 0; i < keys.size(); ++i) param_dict.emplace(keys[i], vals[i]);
    m.updateParams(param_dict);

    ostringstream rowstream;
    rowstream << std::fixed << setprecision(2);
    for (double val: vals) rowstream << val << ",";
    rowstream << std::scientific << setprecision(5);

    if (ord == ORDER::TREE) {
        for (string wc: wc_names) rowstream << eval_wc(m, wc, mubarsq, 0.0) << ",";
    } else if (ord == ORDER::FULL) {
        for (string wc: wc_names) rowstream << eval_wc(m, wc, mubarsq, hb) << ",";
    } else {
        for (string wc: wc_names) rowstream << "[" << eval_wc(m, wc, mubarsq, 0.0) << "," << eval_wc(m, wc, mubarsq, hb) - eval_wc(m, wc, mubarsq, 0.0) << "],";
    }

    string row = rowstream.str();
    return row.substr(0, row.length()-1);
}

void write_to_csv(string fname, Model& m, map<string, vector<double> > p_range_dict, vector<string> wc_names, double mubarsq, ORDER ord) {
    vector<string> keys;
    for(auto it = p_range_dict.begin(); it != p_range_dict.end(); ++it) keys.emplace_back(it->first);
    vector<vector<double> > p_combs = create_param_combs(p_range_dict);

    if (std::filesystem::exists(fname)) {
        std::cout << "Found existing file at path: " << fname << "\n";
        std::cout << "The contents of this file will be overwritten.\n";
        std::filesystem::remove(fname);
    }

    ofstream f1;
    f1.open(fname, ios::out | ios::app);
    ostringstream h_stream;
    for (string key: keys) h_stream << key << ",";
    for (string wc: wc_names) h_stream << wc << ",";

    string header = h_stream.str();
    f1 << header.substr(0, header.length()-1) << "\n";
    for (auto p_comb: p_combs) f1 << create_row(m, keys, p_comb, wc_names, mubarsq, ord) << "\n";
    f1.close();
}
