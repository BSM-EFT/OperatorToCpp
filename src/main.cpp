/**
 * @file write_to_files.cpp
 * @author Suraj Prakash
 * @date 2025-06-02
 * @brief Interface for reading parameter values and WC names from files and
 *        writing WC values to a .csv file
 */

#include "FileIO.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <map>

using std::vector;
using std::unordered_map;
using std::map;
using std::string;

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
