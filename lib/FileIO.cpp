/**
 * @file FileIO.cpp
 * @author Suraj Prakash
 * @date 2025-10-07
 * @brief A suite of utility functions to aid in reading input from and writing output to files
 */

#include "FileIO.h"
#include "OperatorImport.h"

#include <vector>
#include <string>
#include <cctype>
#include <unordered_map>
#include <map>
#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <filesystem>

#define step 0.1

using std::vector;
using std::unordered_map;
using std::map;
using std::string;
using std::isspace;
using std::ifstream;

string trim_right(string s) {
    int i = s.length() - 1;
    while(isspace(s[i])) i--;
    s.erase(s.begin()+i+1,s.end());
    return s;
}

string trim_left(string s) {
    int i = 0;
    while(isspace(s[i])) i++;
    s.erase(s.begin(),s.begin()+i);
    return s;
}

string trim(string s) {
    return trim_left(trim_right(s));
}

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
            p_dict.emplace(trim(p_name), val);
        } else {
            int b_close = rest.find("]");
            int comma = rest.find(",");
            start = std::stod(rest.substr(b_open + 1, comma - b_open - 1));
            end = std::stod(rest.substr(comma + 2, b_close - b_open - 2));
            vector<double> val_vec = create_range(start, end, step);
            p_range_dict.emplace(trim(p_name), val_vec);
        }
    }
}

void read_params(string fname, unordered_map<string, double>& p_dict) {
    ifstream params_file(fname);
    string line, p_name, rest;
    double val;

    while (std::getline(params_file, line)) {
        int sep = line.find(":");
        p_name = line.substr(0,sep);
        rest = line.substr(sep+1,line.length() - p_name.length());
        val = std::stod(rest.substr(1, rest.length()));
        p_dict.emplace(trim(p_name), val);
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
    while (std::getline(wcs_file, line)) wc_names.emplace_back(trim(line)); // trim leading and trailing whitespaces in WC names

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

void check_file_exists(string f) {
    if (std::filesystem::exists(f)) {
        std::cout << "Found existing file at path: " << f << "\n";
        std::cout << "The contents of this file will be overwritten.\n";
        std::filesystem::remove(f);
    }
}
