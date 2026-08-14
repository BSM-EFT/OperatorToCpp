/**
 * @file main.cpp
 * @author Suraj Prakash
 * @date 2026-08-09
 * @brief Example C++ program that creates an instance of the SingletScalarExtension class and evaluates Wilson coefficients
 */

#include "SingletScalarExtension.h"
#include <vector>
#include <string>
#include <complex>
#include <unordered_map>
#include <iostream>

int main() {
    // create a parameter dictionary
    std::unordered_map<std::string, std::complex<double>> param_dict;
    param_dict.emplace("gY", 0.36);
    param_dict.emplace("gL", 0.63);
    param_dict.emplace("lmbd", 0.085);

    param_dict.emplace("M", 1.5);
    param_dict.emplace("kappa", 0.2);
    param_dict.emplace("A", 0.1);
    param_dict.emplace("mu", 0.65);
    param_dict.emplace("lmbdPhi", 0.3);

    // Yukawa couplings
    param_dict.emplace("Yu11", 7e-6);
    param_dict.emplace("Yu22", 3.3e-3);
    param_dict.emplace("Yu33", 0.86);
    
    param_dict.emplace("Yd11", 1.5e-5);
    param_dict.emplace("Yd22", 3e-4);
    param_dict.emplace("Yd33", 0.015);
    
    param_dict.emplace("Ye11", 2.9e-6);
    param_dict.emplace("Ye22", 6e-4);
    param_dict.emplace("Ye33", 0.01);
    
     // Heavy Neutrino Yukawa couplings
    param_dict.emplace("Yn11", 1);
    param_dict.emplace("Yn22", 2);
    param_dict.emplace("Yn33", 3);
    
    // Heavy Neutrino masses
    param_dict.emplace("MNR1", 1e8);
    param_dict.emplace("MNR2", 2e8);
    param_dict.emplace("MNR3", 3e8);

    SingletScalarExtension model1(param_dict, 1.0, true); 

    // compute and print Wilson coefficients values
    std::cout << "cH: {Real = " << model1.cH().real()
              << ", Imag = " << model1.cH().imag() << "}\n";
    std::cout << "cuB: {Real = " << model1.cuB(2,2).real()
              << ", Imag = " << model1.cuB(2,2).imag() << "}\n";

    return 0;
}
