/**
 * @file write_to_files.cpp
 * @author Suraj Prakash
 * @date 2025-07-09
 * @brief Interface for reading parameter values and WC names from files and
 *        writing WC values to a .csv file
 */

#include "FileIO.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <cmath>
#include <iostream>

using std::vector;
using std::unordered_map;
using std::map;
using std::string;
using std::cout;

int main() {

    // create an instance of the model
    Model sb_model;

    // dictionaries to store fixed and variable value parameters
    unordered_map<string, double> par_dict;
    map<string, vector<double> > par_range_dict;

    // read parameter-values and wc-names as input from files
    read_params("./input/params.yaml", par_dict, par_range_dict);
    vector<string> wcs = read_wc_names("./input/coeffs.txt");

    cout << "Input files read without error!" << "\n\n";

    // update the parameters of the model based on the fixed-valued inputs
    sb_model.updateParams(par_dict);
    double mubarsq = pow(par_dict["scale"],2);

    // generate results and store them in a data.csv file
    write_to_csv("./plots/data.csv", sb_model, par_range_dict, wcs, mubarsq);

    cout << "Output written to data.csv!" << "\n";

    return 0;
}
