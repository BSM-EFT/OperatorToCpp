/**
 * @file main.cpp
 * @author Suraj Prakash
 * @date 2026-06-05
 * @brief Example C++ program that creates an instance of the MSSM class and evaluates Wilson coefficients
 */

#include "MSSM.h"
#include <vector>
#include <string>
#include <complex>
#include <unordered_map>
#include <iostream>

int main() {
    // create a parameter dictionary
    std::unordered_map<std::string, std::complex<double>> param_dict;
    param_dict.emplace("g1", 0.37);
    param_dict.emplace("g3", 1.1);
    param_dict.emplace("cgamma",0.01);

    // Bino mass and (right-handed) stop mass (in units of TeV)
    param_dict.emplace("m1", 1.200);
    param_dict.emplace("mut3", 2.000);

    // Yukawa couplings
    param_dict.emplace("Yu11", 0.00001);
    param_dict.emplace("Yu22", 0.007);
    param_dict.emplace("Yu33", 0.9);

    // We set the masses of all other super-partners to be very large, unspecified parmeters remain zero.
    std::vector<std::string> heavy_masses = {
        "m3", "m2", "mPhi", "muTilde",
        "met1", "met2", "met3",
        "mlt1", "mlt2", "mlt3",
        "mqt1", "mqt2", "mqt3",
        "mdt1", "mdt2", "mdt3",
        "mut1", "mut2"
    };

    double i = 0.0;
    for(std::string mass: heavy_masses) {
        param_dict.emplace(mass, 1'000'000.0 + i);
        i += 1000;
    }

    // initialize an instance of the MSSM model with the parameter dictionary,
    // renormalization scale set to 1 TeV and loop contrbutions turned on
    MSSM sb_model(param_dict, 1.0, true); 

    // compute and print Wilson coefficients values
    std::cout << "cG: {Real = " << sb_model.cG().real()
              << ", Imag = " << sb_model.cG().imag() << "}\n";
    std::cout << "cuB: {Real = " << sb_model.cuB(2,2).real()
              << ", Imag = " << sb_model.cuB(2,2).imag() << "}\n";

    return 0;
}
