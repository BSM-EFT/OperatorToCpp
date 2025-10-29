#include "pch.h"
#include "MSSM.h"

MSSM::MSSM(std::unordered_map<std::string, std::complex<double>> params) {
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
    if (params.find("yd11") != params.end()) {
        this->yd[0][0] = params["yd11"];
        this->ydc[0][0] = std::conj(params["yd11"]);
    }
    if (params.find("yd12") != params.end()) {
        this->yd[0][1] = params["yd12"];
        this->ydc[1][0] = std::conj(params["yd12"]);
    }
    if (params.find("yd13") != params.end()) {
        this->yd[0][2] = params["yd13"];
        this->ydc[2][0] = std::conj(params["yd13"]);
    }
    if (params.find("yd21") != params.end()) {
        this->yd[1][0] = params["yd21"];
        this->ydc[0][1] = std::conj(params["yd21"]);
    }
    if (params.find("yd22") != params.end()) {
        this->yd[1][1] = params["yd22"];
        this->ydc[1][1] = std::conj(params["yd22"]);
    }
    if (params.find("yd23") != params.end()) {
        this->yd[1][2] = params["yd23"];
        this->ydc[2][1] = std::conj(params["yd23"]);
    }
    if (params.find("yd31") != params.end()) {
        this->yd[2][0] = params["yd31"];
        this->ydc[0][2] = std::conj(params["yd31"]);
    }
    if (params.find("yd32") != params.end()) {
        this->yd[2][1] = params["yd32"];
        this->ydc[1][2] = std::conj(params["yd32"]);
    }
    if (params.find("yd33") != params.end()) {
        this->yd[2][2] = params["yd33"];
        this->ydc[2][2] = std::conj(params["yd33"]);
    }
    if (params.find("ye11") != params.end()) {
        this->ye[0][0] = params["ye11"];
        this->yec[0][0] = std::conj(params["ye11"]);
    }
    if (params.find("ye12") != params.end()) {
        this->ye[0][1] = params["ye12"];
        this->yec[1][0] = std::conj(params["ye12"]);
    }
    if (params.find("ye13") != params.end()) {
        this->ye[0][2] = params["ye13"];
        this->yec[2][0] = std::conj(params["ye13"]);
    }
    if (params.find("ye21") != params.end()) {
        this->ye[1][0] = params["ye21"];
        this->yec[0][1] = std::conj(params["ye21"]);
    }
    if (params.find("ye22") != params.end()) {
        this->ye[1][1] = params["ye22"];
        this->yec[1][1] = std::conj(params["ye22"]);
    }
    if (params.find("ye23") != params.end()) {
        this->ye[1][2] = params["ye23"];
        this->yec[2][1] = std::conj(params["ye23"]);
    }
    if (params.find("ye31") != params.end()) {
        this->ye[2][0] = params["ye31"];
        this->yec[0][2] = std::conj(params["ye31"]);
    }
    if (params.find("ye32") != params.end()) {
        this->ye[2][1] = params["ye32"];
        this->yec[1][2] = std::conj(params["ye32"]);
    }
    if (params.find("ye33") != params.end()) {
        this->ye[2][2] = params["ye33"];
        this->yec[2][2] = std::conj(params["ye33"]);
    }
    if (params.find("yu11") != params.end()) {
        this->yu[0][0] = params["yu11"];
        this->yuc[0][0] = std::conj(params["yu11"]);
    }
    if (params.find("yu12") != params.end()) {
        this->yu[0][1] = params["yu12"];
        this->yuc[1][0] = std::conj(params["yu12"]);
    }
    if (params.find("yu13") != params.end()) {
        this->yu[0][2] = params["yu13"];
        this->yuc[2][0] = std::conj(params["yu13"]);
    }
    if (params.find("yu21") != params.end()) {
        this->yu[1][0] = params["yu21"];
        this->yuc[0][1] = std::conj(params["yu21"]);
    }
    if (params.find("yu22") != params.end()) {
        this->yu[1][1] = params["yu22"];
        this->yuc[1][1] = std::conj(params["yu22"]);
    }
    if (params.find("yu23") != params.end()) {
        this->yu[1][2] = params["yu23"];
        this->yuc[2][1] = std::conj(params["yu23"]);
    }
    if (params.find("yu31") != params.end()) {
        this->yu[2][0] = params["yu31"];
        this->yuc[0][2] = std::conj(params["yu31"]);
    }
    if (params.find("yu32") != params.end()) {
        this->yu[2][1] = params["yu32"];
        this->yuc[1][2] = std::conj(params["yu32"]);
    }
    if (params.find("yu33") != params.end()) {
        this->yu[2][2] = params["yu33"];
        this->yuc[2][2] = std::conj(params["yu33"]);
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
    if (params.find("yd11") != params.end()) {
        this->yd[0][0] = params["yd11"];
        this->ydc[0][0] = std::conj(params["yd11"]);
    }
    if (params.find("yd12") != params.end()) {
        this->yd[0][1] = params["yd12"];
        this->ydc[1][0] = std::conj(params["yd12"]);
    }
    if (params.find("yd13") != params.end()) {
        this->yd[0][2] = params["yd13"];
        this->ydc[2][0] = std::conj(params["yd13"]);
    }
    if (params.find("yd21") != params.end()) {
        this->yd[1][0] = params["yd21"];
        this->ydc[0][1] = std::conj(params["yd21"]);
    }
    if (params.find("yd22") != params.end()) {
        this->yd[1][1] = params["yd22"];
        this->ydc[1][1] = std::conj(params["yd22"]);
    }
    if (params.find("yd23") != params.end()) {
        this->yd[1][2] = params["yd23"];
        this->ydc[2][1] = std::conj(params["yd23"]);
    }
    if (params.find("yd31") != params.end()) {
        this->yd[2][0] = params["yd31"];
        this->ydc[0][2] = std::conj(params["yd31"]);
    }
    if (params.find("yd32") != params.end()) {
        this->yd[2][1] = params["yd32"];
        this->ydc[1][2] = std::conj(params["yd32"]);
    }
    if (params.find("yd33") != params.end()) {
        this->yd[2][2] = params["yd33"];
        this->ydc[2][2] = std::conj(params["yd33"]);
    }
    if (params.find("ye11") != params.end()) {
        this->ye[0][0] = params["ye11"];
        this->yec[0][0] = std::conj(params["ye11"]);
    }
    if (params.find("ye12") != params.end()) {
        this->ye[0][1] = params["ye12"];
        this->yec[1][0] = std::conj(params["ye12"]);
    }
    if (params.find("ye13") != params.end()) {
        this->ye[0][2] = params["ye13"];
        this->yec[2][0] = std::conj(params["ye13"]);
    }
    if (params.find("ye21") != params.end()) {
        this->ye[1][0] = params["ye21"];
        this->yec[0][1] = std::conj(params["ye21"]);
    }
    if (params.find("ye22") != params.end()) {
        this->ye[1][1] = params["ye22"];
        this->yec[1][1] = std::conj(params["ye22"]);
    }
    if (params.find("ye23") != params.end()) {
        this->ye[1][2] = params["ye23"];
        this->yec[2][1] = std::conj(params["ye23"]);
    }
    if (params.find("ye31") != params.end()) {
        this->ye[2][0] = params["ye31"];
        this->yec[0][2] = std::conj(params["ye31"]);
    }
    if (params.find("ye32") != params.end()) {
        this->ye[2][1] = params["ye32"];
        this->yec[1][2] = std::conj(params["ye32"]);
    }
    if (params.find("ye33") != params.end()) {
        this->ye[2][2] = params["ye33"];
        this->yec[2][2] = std::conj(params["ye33"]);
    }
    if (params.find("yu11") != params.end()) {
        this->yu[0][0] = params["yu11"];
        this->yuc[0][0] = std::conj(params["yu11"]);
    }
    if (params.find("yu12") != params.end()) {
        this->yu[0][1] = params["yu12"];
        this->yuc[1][0] = std::conj(params["yu12"]);
    }
    if (params.find("yu13") != params.end()) {
        this->yu[0][2] = params["yu13"];
        this->yuc[2][0] = std::conj(params["yu13"]);
    }
    if (params.find("yu21") != params.end()) {
        this->yu[1][0] = params["yu21"];
        this->yuc[0][1] = std::conj(params["yu21"]);
    }
    if (params.find("yu22") != params.end()) {
        this->yu[1][1] = params["yu22"];
        this->yuc[1][1] = std::conj(params["yu22"]);
    }
    if (params.find("yu23") != params.end()) {
        this->yu[1][2] = params["yu23"];
        this->yuc[2][1] = std::conj(params["yu23"]);
    }
    if (params.find("yu31") != params.end()) {
        this->yu[2][0] = params["yu31"];
        this->yuc[0][2] = std::conj(params["yu31"]);
    }
    if (params.find("yu32") != params.end()) {
        this->yu[2][1] = params["yu32"];
        this->yuc[1][2] = std::conj(params["yu32"]);
    }
    if (params.find("yu33") != params.end()) {
        this->yu[2][2] = params["yu33"];
        this->yuc[2][2] = std::conj(params["yu33"]);
    }
}
