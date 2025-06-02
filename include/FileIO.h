/**
 * @file FileIO.h
 * @author Suraj Prakash
 * @date 2025-06-02
 * @brief header file correspondig to FileIO.cpp
 */

#include "MSSM.h"

using Model = MSSM;

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <utility>

std::vector<double> create_range(double start, double end, double delta);

void read_params(std::string fname, std::unordered_map<std::string, double>& p_dict, std::map<std::string, std::vector<double> >& p_range_dict);

std::vector<std::vector<double> > create_param_combs(std::map<std::string, std::vector<double> >& p_range_dict);

std::vector<std::string> read_wc_names(std::string fname);

std::pair<std::string,std::vector<int>> split_name_idx(std::string full_name);

double eval_wc(Model m, std::string s, double mubarsq);

std::string create_row(Model& m, std::vector<std::string>& keys, std::vector<double>& vals, std::vector<std::string>& wc_names);

void write_to_csv(Model& m, std::map<std::string, std::vector<double> > p_range_dict, std::vector<std::string> wc_names);
