/**
 * @file write_to_files.cpp
 * @author Suraj Prakash
 * @date 2025-07-11
 * @brief Interface for reading parameter values and WC names from files and
 *        writing WC values to a .csv file
 */

#include "FileIO.h"
#include <cstdlib>
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

    // read parameter-values and wc-names as input from files
    read_params("./input/params-2.yaml", par_dict);
    vector<string> wcs = read_wc_names("./input/coeffs-2.txt");

    cout << "Input files read without error!" << "\n\n";

    // update the parameters of the model based on the fixed-valued inputs
    sb_model.updateParams(par_dict);
    double mubarsq = pow(par_dict["scale"],2);

    vector<vector<double> > par_combs;
    // populate par_combs using random values
    for (int i = 0; i < 50; i++)  {
        vector<double> v;
        double muTilde = (1000 + rand() % 2000)/2000.0;
        int upper = (int)(10-muTilde)*2000;
        double mqt3 = (muTilde*6000 + rand() % (upper*3))/6000.0;
        double mut3 = (muTilde*4000 + rand() % (upper*2))/4000.0;

        v.emplace_back(muTilde);
        v.emplace_back(mqt3);
        v.emplace_back(mut3);
        par_combs.emplace_back(v);
    }


    // generate results and store them in a data.csv file
    write_to_csv(
        "./plots/data-2.csv",
        sb_model,
        par_combs,
        {"muTilde", "mqt3", "mut3"},
        wcs,
        mubarsq,
        ORDER::FULL
    );

    cout << "Output written to data-2.csv!" << "\n";

    return 0;
}
