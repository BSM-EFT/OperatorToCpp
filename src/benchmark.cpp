/**
 * @file benchmark.cpp
 * @author Suraj Prakash
 * @date 2025-10-14
 * @brief Code for estimating time performance of the Wilson coefficient functions in MSSM.h and MSSM.cpp
 */

#include "MSSM.h"
#include <utility>
#include <vector>
#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <ratio>
#include <cmath>

using std::vector;
using std::unordered_map;
using std::string;

int main() {

    MSSM sb_model;

    double yDM = 0.10;
    double gs = 1.1;
    double ytSM = 0.9;
    double mubarsq = 1.05*1.05;

    std::unordered_map<string, double> param_dict;

    param_dict.emplace("g1", yDM*3/(2*sqrt(2)));
    param_dict.emplace("g3", gs);
    param_dict.emplace("cgamma",0.01); // cos(\beta) should not be 0 or 1.

    // Bino mass (in units of TeV), this is M_1 in Eq.(6.1)
    param_dict.emplace("m1", 1.200);

    // (right-handed) stop mass (in units of TeV), this is m_\tilde{t} in Eq.(6.1)
    param_dict.emplace("mut3", 2.000);

    // Higgs mass set to physical mass (units of TeV) (this does not factor into the results)
    param_dict.emplace("mHsq", 0.125*0.125);

    // SM Yukawas
    param_dict.emplace("yu11", 0.00001);
    param_dict.emplace("yu22", 0.007);
    param_dict.emplace("yu33", ytSM);

    // We set the masses of all other superpartners to be very large, unspecified parmeters remain zero.
    vector<string> heavy_masses = {
        "m3", "m2", "mPhi", "muTilde",
        "met1", "met2", "met3",
        "mlt1", "mlt2", "mlt3",
        "mqt1", "mqt2", "mqt3",
        "mdt1", "mdt2", "mdt3",
        "mut1", "mut2"
    };

    double i = 0.0;
    for(string mass: heavy_masses) {
        param_dict.emplace(mass, 1'000'000.0 + i);
        i += 1000;
    }

    sb_model.updateParams(param_dict);

    // evaluation time for given function
    auto eval = [](auto func) {
        const auto t1 = std::chrono::high_resolution_clock::now();
        const auto [name, result] = func();
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << name << " = " << result << ", execution time = " << ms.count() << " ms.\n";
    };

    const auto ti = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < 10; i++) {
        eval([&](){
            return std::pair{"cH", sb_model.cH(mubarsq,0.006332574)};
        });
    }

    const auto tf = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = tf - ti;
    std::cout << "Execution time for 10 calls to cH with parallel EinsSum() and no caching: " << ms.count() << " ms.\n\n";

    return 0;
}
