#include "pch.h"
#include "MSSM.h"
#include <omp.h>

MSSM::MSSM(std::unordered_map<std::string, std::complex<double>> params, double scale, bool loop) {
    this->mubarsq = scale * scale;
    if (!loop) {
        this->hbar = 0.0;
    } else {
        this->hbar = 1/(16 * pow(pi,2));
    }

    if (params.find("cgamma") != params.end()) {
        this->cgamma = params["cgamma"];
    }
    if (params.find("g1") != params.end()) {
        this->g1 = params["g1"];
    }
    if (params.find("g2") != params.end()) {
        this->g2 = params["g2"];
    }
    if (params.find("g3") != params.end()) {
        this->g3 = params["g3"];
    }
    if (params.find("m1") != params.end()) {
        this->m1 = params["m1"];
    }
    if (params.find("m2") != params.end()) {
        this->m2 = params["m2"];
    }
    if (params.find("m3") != params.end()) {
        this->m3 = params["m3"];
    }
    if (params.find("mPhi") != params.end()) {
        this->mPhi = params["mPhi"];
    }
    if (params.find("muTilde") != params.end()) {
        this->muTilde = params["muTilde"];
    }
    if (params.find("mdt1") != params.end()) {
        this->mdt[0] = params["mdt1"];
    }
    if (params.find("mdt2") != params.end()) {
        this->mdt[1] = params["mdt2"];
    }
    if (params.find("mdt3") != params.end()) {
        this->mdt[2] = params["mdt3"];
    }
    if (params.find("met1") != params.end()) {
        this->met[0] = params["met1"];
    }
    if (params.find("met2") != params.end()) {
        this->met[1] = params["met2"];
    }
    if (params.find("met3") != params.end()) {
        this->met[2] = params["met3"];
    }
    if (params.find("mlt1") != params.end()) {
        this->mlt[0] = params["mlt1"];
    }
    if (params.find("mlt2") != params.end()) {
        this->mlt[1] = params["mlt2"];
    }
    if (params.find("mlt3") != params.end()) {
        this->mlt[2] = params["mlt3"];
    }
    if (params.find("mqt1") != params.end()) {
        this->mqt[0] = params["mqt1"];
    }
    if (params.find("mqt2") != params.end()) {
        this->mqt[1] = params["mqt2"];
    }
    if (params.find("mqt3") != params.end()) {
        this->mqt[2] = params["mqt3"];
    }
    if (params.find("mut1") != params.end()) {
        this->mut[0] = params["mut1"];
    }
    if (params.find("mut2") != params.end()) {
        this->mut[1] = params["mut2"];
    }
    if (params.find("mut3") != params.end()) {
        this->mut[2] = params["mut3"];
    }
    if (params.find("ad11") != params.end()) {
        this->ad[0][0] = params["ad11"];
        this->adc[0][0] = std::conj(params["ad11"]);
    }
    if (params.find("ad12") != params.end()) {
        this->ad[0][1] = params["ad12"];
        this->adc[1][0] = std::conj(params["ad12"]);
    }
    if (params.find("ad13") != params.end()) {
        this->ad[0][2] = params["ad13"];
        this->adc[2][0] = std::conj(params["ad13"]);
    }
    if (params.find("ad21") != params.end()) {
        this->ad[1][0] = params["ad21"];
        this->adc[0][1] = std::conj(params["ad21"]);
    }
    if (params.find("ad22") != params.end()) {
        this->ad[1][1] = params["ad22"];
        this->adc[1][1] = std::conj(params["ad22"]);
    }
    if (params.find("ad23") != params.end()) {
        this->ad[1][2] = params["ad23"];
        this->adc[2][1] = std::conj(params["ad23"]);
    }
    if (params.find("ad31") != params.end()) {
        this->ad[2][0] = params["ad31"];
        this->adc[0][2] = std::conj(params["ad31"]);
    }
    if (params.find("ad32") != params.end()) {
        this->ad[2][1] = params["ad32"];
        this->adc[1][2] = std::conj(params["ad32"]);
    }
    if (params.find("ad33") != params.end()) {
        this->ad[2][2] = params["ad33"];
        this->adc[2][2] = std::conj(params["ad33"]);
    }
    if (params.find("ae11") != params.end()) {
        this->ae[0][0] = params["ae11"];
        this->aec[0][0] = std::conj(params["ae11"]);
    }
    if (params.find("ae12") != params.end()) {
        this->ae[0][1] = params["ae12"];
        this->aec[1][0] = std::conj(params["ae12"]);
    }
    if (params.find("ae13") != params.end()) {
        this->ae[0][2] = params["ae13"];
        this->aec[2][0] = std::conj(params["ae13"]);
    }
    if (params.find("ae21") != params.end()) {
        this->ae[1][0] = params["ae21"];
        this->aec[0][1] = std::conj(params["ae21"]);
    }
    if (params.find("ae22") != params.end()) {
        this->ae[1][1] = params["ae22"];
        this->aec[1][1] = std::conj(params["ae22"]);
    }
    if (params.find("ae23") != params.end()) {
        this->ae[1][2] = params["ae23"];
        this->aec[2][1] = std::conj(params["ae23"]);
    }
    if (params.find("ae31") != params.end()) {
        this->ae[2][0] = params["ae31"];
        this->aec[0][2] = std::conj(params["ae31"]);
    }
    if (params.find("ae32") != params.end()) {
        this->ae[2][1] = params["ae32"];
        this->aec[1][2] = std::conj(params["ae32"]);
    }
    if (params.find("ae33") != params.end()) {
        this->ae[2][2] = params["ae33"];
        this->aec[2][2] = std::conj(params["ae33"]);
    }
    if (params.find("au11") != params.end()) {
        this->au[0][0] = params["au11"];
        this->auc[0][0] = std::conj(params["au11"]);
    }
    if (params.find("au12") != params.end()) {
        this->au[0][1] = params["au12"];
        this->auc[1][0] = std::conj(params["au12"]);
    }
    if (params.find("au13") != params.end()) {
        this->au[0][2] = params["au13"];
        this->auc[2][0] = std::conj(params["au13"]);
    }
    if (params.find("au21") != params.end()) {
        this->au[1][0] = params["au21"];
        this->auc[0][1] = std::conj(params["au21"]);
    }
    if (params.find("au22") != params.end()) {
        this->au[1][1] = params["au22"];
        this->auc[1][1] = std::conj(params["au22"]);
    }
    if (params.find("au23") != params.end()) {
        this->au[1][2] = params["au23"];
        this->auc[2][1] = std::conj(params["au23"]);
    }
    if (params.find("au31") != params.end()) {
        this->au[2][0] = params["au31"];
        this->auc[0][2] = std::conj(params["au31"]);
    }
    if (params.find("au32") != params.end()) {
        this->au[2][1] = params["au32"];
        this->auc[1][2] = std::conj(params["au32"]);
    }
    if (params.find("au33") != params.end()) {
        this->au[2][2] = params["au33"];
        this->auc[2][2] = std::conj(params["au33"]);
    }
    if (params.find("Yd11") != params.end()) {
        this->Yd[0][0] = params["Yd11"];
        this->Ydc[0][0] = std::conj(params["Yd11"]);
    }
    if (params.find("Yd12") != params.end()) {
        this->Yd[0][1] = params["Yd12"];
        this->Ydc[1][0] = std::conj(params["Yd12"]);
    }
    if (params.find("Yd13") != params.end()) {
        this->Yd[0][2] = params["Yd13"];
        this->Ydc[2][0] = std::conj(params["Yd13"]);
    }
    if (params.find("Yd21") != params.end()) {
        this->Yd[1][0] = params["Yd21"];
        this->Ydc[0][1] = std::conj(params["Yd21"]);
    }
    if (params.find("Yd22") != params.end()) {
        this->Yd[1][1] = params["Yd22"];
        this->Ydc[1][1] = std::conj(params["Yd22"]);
    }
    if (params.find("Yd23") != params.end()) {
        this->Yd[1][2] = params["Yd23"];
        this->Ydc[2][1] = std::conj(params["Yd23"]);
    }
    if (params.find("Yd31") != params.end()) {
        this->Yd[2][0] = params["Yd31"];
        this->Ydc[0][2] = std::conj(params["Yd31"]);
    }
    if (params.find("Yd32") != params.end()) {
        this->Yd[2][1] = params["Yd32"];
        this->Ydc[1][2] = std::conj(params["Yd32"]);
    }
    if (params.find("Yd33") != params.end()) {
        this->Yd[2][2] = params["Yd33"];
        this->Ydc[2][2] = std::conj(params["Yd33"]);
    }
    if (params.find("Ye11") != params.end()) {
        this->Ye[0][0] = params["Ye11"];
        this->Yec[0][0] = std::conj(params["Ye11"]);
    }
    if (params.find("Ye12") != params.end()) {
        this->Ye[0][1] = params["Ye12"];
        this->Yec[1][0] = std::conj(params["Ye12"]);
    }
    if (params.find("Ye13") != params.end()) {
        this->Ye[0][2] = params["Ye13"];
        this->Yec[2][0] = std::conj(params["Ye13"]);
    }
    if (params.find("Ye21") != params.end()) {
        this->Ye[1][0] = params["Ye21"];
        this->Yec[0][1] = std::conj(params["Ye21"]);
    }
    if (params.find("Ye22") != params.end()) {
        this->Ye[1][1] = params["Ye22"];
        this->Yec[1][1] = std::conj(params["Ye22"]);
    }
    if (params.find("Ye23") != params.end()) {
        this->Ye[1][2] = params["Ye23"];
        this->Yec[2][1] = std::conj(params["Ye23"]);
    }
    if (params.find("Ye31") != params.end()) {
        this->Ye[2][0] = params["Ye31"];
        this->Yec[0][2] = std::conj(params["Ye31"]);
    }
    if (params.find("Ye32") != params.end()) {
        this->Ye[2][1] = params["Ye32"];
        this->Yec[1][2] = std::conj(params["Ye32"]);
    }
    if (params.find("Ye33") != params.end()) {
        this->Ye[2][2] = params["Ye33"];
        this->Yec[2][2] = std::conj(params["Ye33"]);
    }
    if (params.find("Yu11") != params.end()) {
        this->Yu[0][0] = params["Yu11"];
        this->Yuc[0][0] = std::conj(params["Yu11"]);
    }
    if (params.find("Yu12") != params.end()) {
        this->Yu[0][1] = params["Yu12"];
        this->Yuc[1][0] = std::conj(params["Yu12"]);
    }
    if (params.find("Yu13") != params.end()) {
        this->Yu[0][2] = params["Yu13"];
        this->Yuc[2][0] = std::conj(params["Yu13"]);
    }
    if (params.find("Yu21") != params.end()) {
        this->Yu[1][0] = params["Yu21"];
        this->Yuc[0][1] = std::conj(params["Yu21"]);
    }
    if (params.find("Yu22") != params.end()) {
        this->Yu[1][1] = params["Yu22"];
        this->Yuc[1][1] = std::conj(params["Yu22"]);
    }
    if (params.find("Yu23") != params.end()) {
        this->Yu[1][2] = params["Yu23"];
        this->Yuc[2][1] = std::conj(params["Yu23"]);
    }
    if (params.find("Yu31") != params.end()) {
        this->Yu[2][0] = params["Yu31"];
        this->Yuc[0][2] = std::conj(params["Yu31"]);
    }
    if (params.find("Yu32") != params.end()) {
        this->Yu[2][1] = params["Yu32"];
        this->Yuc[1][2] = std::conj(params["Yu32"]);
    }
    if (params.find("Yu33") != params.end()) {
        this->Yu[2][2] = params["Yu33"];
        this->Yuc[2][2] = std::conj(params["Yu33"]);
    }
}

void MSSM::updateParams(std::unordered_map<std::string, std::complex<double>> params) {
    if (params.find("cgamma") != params.end()) {
        this->cgamma = params["cgamma"];
    }
    if (params.find("g1") != params.end()) {
        this->g1 = params["g1"];
    }
    if (params.find("g2") != params.end()) {
        this->g2 = params["g2"];
    }
    if (params.find("g3") != params.end()) {
        this->g3 = params["g3"];
    }
    if (params.find("m1") != params.end()) {
        this->m1 = params["m1"];
    }
    if (params.find("m2") != params.end()) {
        this->m2 = params["m2"];
    }
    if (params.find("m3") != params.end()) {
        this->m3 = params["m3"];
    }
    if (params.find("mPhi") != params.end()) {
        this->mPhi = params["mPhi"];
    }
    if (params.find("muTilde") != params.end()) {
        this->muTilde = params["muTilde"];
    }
    if (params.find("mdt1") != params.end()) {
        this->mdt[0] = params["mdt1"];
    }
    if (params.find("mdt2") != params.end()) {
        this->mdt[1] = params["mdt2"];
    }
    if (params.find("mdt3") != params.end()) {
        this->mdt[2] = params["mdt3"];
    }
    if (params.find("met1") != params.end()) {
        this->met[0] = params["met1"];
    }
    if (params.find("met2") != params.end()) {
        this->met[1] = params["met2"];
    }
    if (params.find("met3") != params.end()) {
        this->met[2] = params["met3"];
    }
    if (params.find("mlt1") != params.end()) {
        this->mlt[0] = params["mlt1"];
    }
    if (params.find("mlt2") != params.end()) {
        this->mlt[1] = params["mlt2"];
    }
    if (params.find("mlt3") != params.end()) {
        this->mlt[2] = params["mlt3"];
    }
    if (params.find("mqt1") != params.end()) {
        this->mqt[0] = params["mqt1"];
    }
    if (params.find("mqt2") != params.end()) {
        this->mqt[1] = params["mqt2"];
    }
    if (params.find("mqt3") != params.end()) {
        this->mqt[2] = params["mqt3"];
    }
    if (params.find("mut1") != params.end()) {
        this->mut[0] = params["mut1"];
    }
    if (params.find("mut2") != params.end()) {
        this->mut[1] = params["mut2"];
    }
    if (params.find("mut3") != params.end()) {
        this->mut[2] = params["mut3"];
    }
    if (params.find("ad11") != params.end()) {
        this->ad[0][0] = params["ad11"];
        this->adc[0][0] = std::conj(params["ad11"]);
    }
    if (params.find("ad12") != params.end()) {
        this->ad[0][1] = params["ad12"];
        this->adc[1][0] = std::conj(params["ad12"]);
    }
    if (params.find("ad13") != params.end()) {
        this->ad[0][2] = params["ad13"];
        this->adc[2][0] = std::conj(params["ad13"]);
    }
    if (params.find("ad21") != params.end()) {
        this->ad[1][0] = params["ad21"];
        this->adc[0][1] = std::conj(params["ad21"]);
    }
    if (params.find("ad22") != params.end()) {
        this->ad[1][1] = params["ad22"];
        this->adc[1][1] = std::conj(params["ad22"]);
    }
    if (params.find("ad23") != params.end()) {
        this->ad[1][2] = params["ad23"];
        this->adc[2][1] = std::conj(params["ad23"]);
    }
    if (params.find("ad31") != params.end()) {
        this->ad[2][0] = params["ad31"];
        this->adc[0][2] = std::conj(params["ad31"]);
    }
    if (params.find("ad32") != params.end()) {
        this->ad[2][1] = params["ad32"];
        this->adc[1][2] = std::conj(params["ad32"]);
    }
    if (params.find("ad33") != params.end()) {
        this->ad[2][2] = params["ad33"];
        this->adc[2][2] = std::conj(params["ad33"]);
    }
    if (params.find("ae11") != params.end()) {
        this->ae[0][0] = params["ae11"];
        this->aec[0][0] = std::conj(params["ae11"]);
    }
    if (params.find("ae12") != params.end()) {
        this->ae[0][1] = params["ae12"];
        this->aec[1][0] = std::conj(params["ae12"]);
    }
    if (params.find("ae13") != params.end()) {
        this->ae[0][2] = params["ae13"];
        this->aec[2][0] = std::conj(params["ae13"]);
    }
    if (params.find("ae21") != params.end()) {
        this->ae[1][0] = params["ae21"];
        this->aec[0][1] = std::conj(params["ae21"]);
    }
    if (params.find("ae22") != params.end()) {
        this->ae[1][1] = params["ae22"];
        this->aec[1][1] = std::conj(params["ae22"]);
    }
    if (params.find("ae23") != params.end()) {
        this->ae[1][2] = params["ae23"];
        this->aec[2][1] = std::conj(params["ae23"]);
    }
    if (params.find("ae31") != params.end()) {
        this->ae[2][0] = params["ae31"];
        this->aec[0][2] = std::conj(params["ae31"]);
    }
    if (params.find("ae32") != params.end()) {
        this->ae[2][1] = params["ae32"];
        this->aec[1][2] = std::conj(params["ae32"]);
    }
    if (params.find("ae33") != params.end()) {
        this->ae[2][2] = params["ae33"];
        this->aec[2][2] = std::conj(params["ae33"]);
    }
    if (params.find("au11") != params.end()) {
        this->au[0][0] = params["au11"];
        this->auc[0][0] = std::conj(params["au11"]);
    }
    if (params.find("au12") != params.end()) {
        this->au[0][1] = params["au12"];
        this->auc[1][0] = std::conj(params["au12"]);
    }
    if (params.find("au13") != params.end()) {
        this->au[0][2] = params["au13"];
        this->auc[2][0] = std::conj(params["au13"]);
    }
    if (params.find("au21") != params.end()) {
        this->au[1][0] = params["au21"];
        this->auc[0][1] = std::conj(params["au21"]);
    }
    if (params.find("au22") != params.end()) {
        this->au[1][1] = params["au22"];
        this->auc[1][1] = std::conj(params["au22"]);
    }
    if (params.find("au23") != params.end()) {
        this->au[1][2] = params["au23"];
        this->auc[2][1] = std::conj(params["au23"]);
    }
    if (params.find("au31") != params.end()) {
        this->au[2][0] = params["au31"];
        this->auc[0][2] = std::conj(params["au31"]);
    }
    if (params.find("au32") != params.end()) {
        this->au[2][1] = params["au32"];
        this->auc[1][2] = std::conj(params["au32"]);
    }
    if (params.find("au33") != params.end()) {
        this->au[2][2] = params["au33"];
        this->auc[2][2] = std::conj(params["au33"]);
    }
    if (params.find("Yd11") != params.end()) {
        this->Yd[0][0] = params["Yd11"];
        this->Ydc[0][0] = std::conj(params["Yd11"]);
    }
    if (params.find("Yd12") != params.end()) {
        this->Yd[0][1] = params["Yd12"];
        this->Ydc[1][0] = std::conj(params["Yd12"]);
    }
    if (params.find("Yd13") != params.end()) {
        this->Yd[0][2] = params["Yd13"];
        this->Ydc[2][0] = std::conj(params["Yd13"]);
    }
    if (params.find("Yd21") != params.end()) {
        this->Yd[1][0] = params["Yd21"];
        this->Ydc[0][1] = std::conj(params["Yd21"]);
    }
    if (params.find("Yd22") != params.end()) {
        this->Yd[1][1] = params["Yd22"];
        this->Ydc[1][1] = std::conj(params["Yd22"]);
    }
    if (params.find("Yd23") != params.end()) {
        this->Yd[1][2] = params["Yd23"];
        this->Ydc[2][1] = std::conj(params["Yd23"]);
    }
    if (params.find("Yd31") != params.end()) {
        this->Yd[2][0] = params["Yd31"];
        this->Ydc[0][2] = std::conj(params["Yd31"]);
    }
    if (params.find("Yd32") != params.end()) {
        this->Yd[2][1] = params["Yd32"];
        this->Ydc[1][2] = std::conj(params["Yd32"]);
    }
    if (params.find("Yd33") != params.end()) {
        this->Yd[2][2] = params["Yd33"];
        this->Ydc[2][2] = std::conj(params["Yd33"]);
    }
    if (params.find("Ye11") != params.end()) {
        this->Ye[0][0] = params["Ye11"];
        this->Yec[0][0] = std::conj(params["Ye11"]);
    }
    if (params.find("Ye12") != params.end()) {
        this->Ye[0][1] = params["Ye12"];
        this->Yec[1][0] = std::conj(params["Ye12"]);
    }
    if (params.find("Ye13") != params.end()) {
        this->Ye[0][2] = params["Ye13"];
        this->Yec[2][0] = std::conj(params["Ye13"]);
    }
    if (params.find("Ye21") != params.end()) {
        this->Ye[1][0] = params["Ye21"];
        this->Yec[0][1] = std::conj(params["Ye21"]);
    }
    if (params.find("Ye22") != params.end()) {
        this->Ye[1][1] = params["Ye22"];
        this->Yec[1][1] = std::conj(params["Ye22"]);
    }
    if (params.find("Ye23") != params.end()) {
        this->Ye[1][2] = params["Ye23"];
        this->Yec[2][1] = std::conj(params["Ye23"]);
    }
    if (params.find("Ye31") != params.end()) {
        this->Ye[2][0] = params["Ye31"];
        this->Yec[0][2] = std::conj(params["Ye31"]);
    }
    if (params.find("Ye32") != params.end()) {
        this->Ye[2][1] = params["Ye32"];
        this->Yec[1][2] = std::conj(params["Ye32"]);
    }
    if (params.find("Ye33") != params.end()) {
        this->Ye[2][2] = params["Ye33"];
        this->Yec[2][2] = std::conj(params["Ye33"]);
    }
    if (params.find("Yu11") != params.end()) {
        this->Yu[0][0] = params["Yu11"];
        this->Yuc[0][0] = std::conj(params["Yu11"]);
    }
    if (params.find("Yu12") != params.end()) {
        this->Yu[0][1] = params["Yu12"];
        this->Yuc[1][0] = std::conj(params["Yu12"]);
    }
    if (params.find("Yu13") != params.end()) {
        this->Yu[0][2] = params["Yu13"];
        this->Yuc[2][0] = std::conj(params["Yu13"]);
    }
    if (params.find("Yu21") != params.end()) {
        this->Yu[1][0] = params["Yu21"];
        this->Yuc[0][1] = std::conj(params["Yu21"]);
    }
    if (params.find("Yu22") != params.end()) {
        this->Yu[1][1] = params["Yu22"];
        this->Yuc[1][1] = std::conj(params["Yu22"]);
    }
    if (params.find("Yu23") != params.end()) {
        this->Yu[1][2] = params["Yu23"];
        this->Yuc[2][1] = std::conj(params["Yu23"]);
    }
    if (params.find("Yu31") != params.end()) {
        this->Yu[2][0] = params["Yu31"];
        this->Yuc[0][2] = std::conj(params["Yu31"]);
    }
    if (params.find("Yu32") != params.end()) {
        this->Yu[2][1] = params["Yu32"];
        this->Yuc[1][2] = std::conj(params["Yu32"]);
    }
    if (params.find("Yu33") != params.end()) {
        this->Yu[2][2] = params["Yu33"];
        this->Yuc[2][2] = std::conj(params["Yu33"]);
    }
}

double MSSM::getScale() {
    return sqrt(this->mubarsq);
}

void MSSM::setScale(double scale) {
    this->mubarsq = scale * scale;
}

void MSSM::loopContributions(bool loop) {
    if (!loop) {
        this->hbar = 0.0;
    } else {
        this->hbar = 1/(16 * pow(pi,2));
    }
}

std::unordered_map<std::string, std::complex<double> > MSSM::getParams() {
    std::unordered_map<std::string, std::complex<double> > param_dict = {
        {"cgamma", this->cgamma},
        {"g1", this->g1},
        {"g2", this->g2},
        {"g3", this->g3},
        {"m1", this->m1},
        {"m2", this->m2},
        {"m3", this->m3},
        {"mPhi", this->mPhi},
        {"muTilde", this->muTilde},
        {"mdt1", this->mdt[0]},
        {"mdt2", this->mdt[1]},
        {"mdt3", this->mdt[2]},
        {"met1", this->met[0]},
        {"met2", this->met[1]},
        {"met3", this->met[2]},
        {"mlt1", this->mlt[0]},
        {"mlt2", this->mlt[1]},
        {"mlt3", this->mlt[2]},
        {"mqt1", this->mqt[0]},
        {"mqt2", this->mqt[1]},
        {"mqt3", this->mqt[2]},
        {"mut1", this->mut[0]},
        {"mut2", this->mut[1]},
        {"mut3", this->mut[2]},
        {"ad11", this->ad[0][0]},
        {"ad12", this->ad[0][1]},
        {"ad13", this->ad[0][2]},
        {"ad22", this->ad[1][1]},
        {"ad23", this->ad[1][2]},
        {"ad33", this->ad[2][2]},
        {"ae11", this->ae[0][0]},
        {"ae12", this->ae[0][1]},
        {"ae13", this->ae[0][2]},
        {"ae22", this->ae[1][1]},
        {"ae23", this->ae[1][2]},
        {"ae33", this->ae[2][2]},
        {"au11", this->au[0][0]},
        {"au12", this->au[0][1]},
        {"au13", this->au[0][2]},
        {"au22", this->au[1][1]},
        {"au23", this->au[1][2]},
        {"au33", this->au[2][2]},
        {"Yd11", this->Yd[0][0]},
        {"Yd12", this->Yd[0][1]},
        {"Yd13", this->Yd[0][2]},
        {"Yd22", this->Yd[1][1]},
        {"Yd23", this->Yd[1][2]},
        {"Yd33", this->Yd[2][2]},
        {"Ye11", this->Ye[0][0]},
        {"Ye12", this->Ye[0][1]},
        {"Ye13", this->Ye[0][2]},
        {"Ye22", this->Ye[1][1]},
        {"Ye23", this->Ye[1][2]},
        {"Ye33", this->Ye[2][2]},
        {"Yu11", this->Yu[0][0]},
        {"Yu12", this->Yu[0][1]},
        {"Yu13", this->Yu[0][2]},
        {"Yu22", this->Yu[1][1]},
        {"Yu23", this->Yu[1][2]},
        {"Yu33", this->Yu[2][2]},
    };
    return param_dict;
}

std::map<std::string, std::complex<double> > MSSM::batch_eval(const std::vector<Task>& tasks) {
    int n = tasks.size();
    std::vector<std::complex<double> > results_temp(n);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        results_temp[i] = tasks[i].work();
    }

    std::map<std::string, std::complex<double> > results;
    for (int i = 0; i < n; ++i) {
        results[tasks[i].name] = results_temp[i];
    }

    return results;
}
