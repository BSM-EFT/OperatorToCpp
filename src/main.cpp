/**
 * @file main.cpp
 * @author Suraj Prakash
 * @date 2025-10-07
 * @brief Interface for reading parameter values and WC names from files and
 *        writing WC values to a .csv file
 */

#include "MSSM.h"
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
    MSSM sb_model;

    // dictionaries to store fixed and variable value parameters
    unordered_map<string, double> par_dict;
    map<string, vector<double> > par_range_dict;

    // read parameter-values and wc-names as input from files
    read_params("./input/params-test.yaml", par_dict, par_range_dict);
    vector<string> wcs = read_wc_names("./input/coeffs-test.txt");

    cout << "Input files read without error!" << "\n\n";

    // update the parameters of the model based on the fixed-valued inputs
    sb_model.updateParams(par_dict);
    double mubarsq = pow(par_dict["scale"],2);

    // generate results and store them in a .yaml file
    write_to_yaml(
        "./plots/data-test.yaml",
        sb_model,
        par_dict,
        {"mPhi","muTilde"},
        wcs,
        mubarsq,
        ORDER::SPLIT
    );

    cout << "Output written to data-test.yaml!" << "\n";

    return 0;
}
