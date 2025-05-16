/**
 * @file write_to_files.cpp
 * @author Suraj Prakash
 * @date 2025-05-16
 * @brief Code for writing (m1, mut3, WC) values to file using MSSM.h and MSSM.cpp
 */

#include "MSSM.h"
#include <vector>
#include <print>
#include <string>
#include <map>
#include <cmath>
#include <format>
#include <iostream>
#include <fstream>

using std::vector;
using std::println;
using std::map;
using std::string;
using std::format;
using std::ofstream;
using std::ios;

int main() {

    MSSM sb_model;

    double g1 = 0.37;
    double gs = 1.1;
    double mubarsq = 1.0*1.0; // We set the renormalization scale close to 1 TeV

    map<string, double> param_dict;

    param_dict.emplace("g1", g1);
    param_dict.emplace("g3", gs);
    param_dict.emplace("cgamma",0.01); // cos(\gamma) should not be 0 or 1.

    // Higgs mass set to physical mass (units of TeV) (this does not factor into the results)
    param_dict.emplace("mHsq", 0.125*0.125);

    // SM Yukawas
    param_dict.emplace("yu11", 0.00001);
    param_dict.emplace("yu22", 0.007);
    param_dict.emplace("yu33", 0.9);

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

    // a lambda that computes WCs for given (mut3, m1) values and creates a string containing all three
    auto compute_and_write = [&sb_model](map<string, double> p_dict, double mut3, double m1, auto func){
        p_dict.emplace("mut3", mut3); // (right-handed) stop mass (in units of TeV)
        p_dict.emplace("m1", m1); // Bino mass (in units of TeV)
        sb_model.updateParams(p_dict);

        auto val = func();
        string line = format("{:.2f} {:.2f} {:.5e}", mut3, m1, val);

        return line;
    };

    // a lambda for writing data to files corresponding to specific WC functions
    auto write_wc = [&param_dict, &compute_and_write](string func_name, auto func){
        string f_name = format("./plots/{}.txt", func_name);

        ofstream file1;
        file1.open(f_name, ios::out | ios::app);
        
        for (double mut3 = 0.3; mut3 < 2.8; mut3 += 0.1) {
            for (double m1 = 0.3; m1 < 2.8; m1 += 0.1) {
                file1 << compute_and_write(param_dict, mut3, m1, func) << "\n";
            }
        }

        file1.close();
    };

    // function calls to write data to files for 2d plots
    write_wc("cG", [&sb_model, &mubarsq](){ return sb_model.cG(mubarsq);});
    write_wc("cuG_33", [&sb_model, &mubarsq](){ return sb_model.cuG(2,2,mubarsq);});
    write_wc("cqu1_1133", [&sb_model, &mubarsq](){ return sb_model.cqu1(0,0,2,2,mubarsq);});
    write_wc("cuu_3333", [&sb_model, &mubarsq](){ return sb_model.cuu(2,2,2,2,mubarsq);});
    write_wc("cqq1_3333", [&sb_model, &mubarsq](){ return sb_model.cqq1(2,2,2,2,mubarsq);});
    write_wc("cqd1_3311", [&sb_model, &mubarsq](){ return sb_model.cqd1(2,2,0,0,mubarsq);});
    write_wc("cqu8_3311", [&sb_model, &mubarsq](){ return sb_model.cqu8(2,2,0,0,mubarsq);});
    write_wc("cqu8_1133", [&sb_model, &mubarsq](){ return sb_model.cqu8(0,0,2,2,mubarsq);});
    

    //  write data to a yaml file to creating bar-chart for multiple benchmark points
    auto create_row = [&sb_model, &mubarsq](map<string, double> p_dict, double mut3, double m1){
        p_dict.emplace("mut3", mut3); // (right-handed) stop mass (in units of TeV)
        p_dict.emplace("m1", m1); // Bino mass (in units of TeV)
        sb_model.updateParams(p_dict);

        string line = format("{:.1f},{:.1f},{:.5e},{:.5e},{:.5e},{:.5e},{:.5e},{:.5e},{:.5e},{:.5e}", 
            mut3, m1, 
            sb_model.cG(mubarsq),
            sb_model.cuG(2,2,mubarsq),
            sb_model.cqu1(0,0,2,2,mubarsq),
            sb_model.cuu(2,2,2,2,mubarsq),
            sb_model.cqq1(2,2,2,2,mubarsq),
            sb_model.cqd1(2,2,0,0,mubarsq),
            sb_model.cqu8(2,2,0,0,mubarsq),
            sb_model.cqu8(0,0,2,2,mubarsq)
        );

        return line;
    };

    string fname = format("./plots/barplot-data.csv");
    ofstream f1;
    f1.open(fname, ios::out | ios::app);
    string first_row = format("mut3,m1,cG,cuG_33,cqu1_1133,cuu_3333,cqq1_3333,cqd1_3311,cqu8_3311,cqu8_1133");
    f1 << first_row << "\n";
    f1 << create_row(param_dict, 2.0, 1.5)  << "\n";    // mut3 = 2.0, m1 = 1.5
    f1 << create_row(param_dict, 1.6, 1.5)  << "\n";    // mut3 = 1.6, m1 = 1.5
    f1 << create_row(param_dict, 1.6, 0.5)  << "\n";    // mut3 = 1.6, m1 = 0.5
    f1.close();

    return 0;
}
