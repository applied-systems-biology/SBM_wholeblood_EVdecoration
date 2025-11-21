/*
 * File:   PI_SBM_Solver.cpp
 *
 * Original work Copyright (c) <2014-2021> by Teresa Lehnert and Maria Prauße
 * Modifications Copyright (c) <2023-2025> by Anastasia Solomatina and Paul Rudolph
 *
 * Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
 * https://www.leibniz-hki.de/en/applied-systems-biology.html
 * HKI-Center for Systems Biology of Infection
 * Leibniz Institute for Natural Product Research and Infection Biology -
 * Hans Knöll Insitute (HKI)
 * Adolf-Reichwein-Straße 23, 07745 Jena, Germany
 *
 * This code is licensed under BSD 2-Clause
 * See the LICENSE file provided with this code for the full license.

 */

#include "model/PI_SBM_Solver.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <set>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <utils/io_utils.h>
#include <vector>
const bool PI_SBM_Solver::REINC_MODE = false;

void PI_SBM_Solver::setSimulationProperties(const sbm::util::SimulationParameters& simulation_parameters) {
    time = 0.0;
    time_step = simulation_parameters.timestep;
    record_time = 0.0;
    record_data_time_step = time_step;
    run_time = simulation_parameters.timesteps;
    no_l = 10; // no_l-1
    no_d = 10; // no_d-1
}

std::vector<double>& PI_SBM_Solver::getStartConcentration() {
    start_concentrations.clear();
    start_concentrations.push_back(temp_Y0);
    start_concentrations.push_back(start_no_M);
    start_concentrations.push_back(start_no_G);
    return start_concentrations;
}

void PI_SBM_Solver::getCollectedStates(std::vector<std::vector<double>>& Y, std::vector<std::vector<std::vector<double>>>& M, std::vector<std::vector<std::vector<double>>>& G, std::vector<std::vector<double>>& R, std::vector<std::vector<double>>& ICk) {
    Y = collected_Y_states_no;
    M = collected_M_states_no;
    G = collected_G_states_no;
    R = collected_no_KE;
    ICk = collected_killedinIC_states_no;
}
void PI_SBM_Solver::set_Parameter(std::vector<double>& parameter) {
    model_parameter.reserve(parameter.size());
    model_parameter.clear();

    for (int i = 0; i < parameter.size(); i++) {
        model_parameter.push_back(parameter[i]);
    }
    setModel();

}

void PI_SBM_Solver::setModel() {
    // set all model variables that must be initialized

    temp_YAE = model_parameter[0]; // C_AE
    temp_YR = model_parameter[1];
    temp_YRA = model_parameter[2];
    temp_YRK = model_parameter[3];
    temp_YKE = model_parameter[4];
    temp_YAG = 0.0;
    temp_YAM = 0.0;
    temp_YKG = 0.0;
    temp_YKM = 0.0;

    start_no_Y_states = temp_YAE + temp_YKE + temp_YAG + temp_YAM; // reactive Y states without resistance and intracellular killed

    if (REINC_MODE == true) {
        temp_Y0 = 1000000.0;
    } else {
        temp_Y0 = (double)temp_YAE + (double)temp_YR + (double)temp_YKE;
    }

    // ini tempM and G
    temp_M_states.clear();
    temp_G_states.clear();
    std::vector<double> temp;
    for (int i = 0; i < no_l; i++) { //no_l --> number of living candida cells
        temp.clear();
        for (int j = 0; j < no_d; j++) { //no_d --> number of dead candida cells
            temp.push_back(0);
        }
        temp_M_states.push_back(temp);
        temp_G_states.push_back(temp);
    }

    temp_M_states[0][0] = start_no_M = model_parameter[5];
    temp_G_states[0][0] = start_no_G = model_parameter[105];

    prelim_M_states = temp_M_states;
    prelim_G_states = temp_G_states;

    prelim_YAE = temp_YAE;
    prelim_YKE = temp_YKE;
    prelim_YR = temp_YR;
    prelim_YRA = temp_YRA;
    prelim_YRK = temp_YRK;
    prelim_YAM = temp_YAM;
    prelim_YAG = temp_YAG;
    prelim_YKM = temp_YKM;
    prelim_YKG = temp_YKG;

    toAddYAE = 0.0;
    toAddYKE = 0.0;
    toAddYAG = 0.0;
    toAddYAM = 0.0;
    toAddYRA = 0.0;
    toAddYRK = 0.0, toAddYKG = 0.0;
    toAddYKM = 0.0;

    occupied_G_states.clear();
    occupied_M_states.clear();
    occupied_GA_states.clear();
    occupied_MA_states.clear();
    int x = 0;
    occupied_G_states.insert(x);
    occupied_M_states.insert(x);

    phiG = model_parameter[model_parameter.size() - 9];
    phiGstar = model_parameter[model_parameter.size() - 8];
    phiM = model_parameter[model_parameter.size() - 7];
    kappaM = model_parameter[model_parameter.size() - 6];
    kappaG = model_parameter[model_parameter.size() - 5];
    rho = model_parameter[model_parameter.size() - 4];
    gamma = model_parameter[model_parameter.size() - 3];
    kappaEG = model_parameter[model_parameter.size() - 2];
    gammaR = model_parameter[model_parameter.size() - 1];

    sum_Individuals = start_no_M + start_no_G + start_no_Y_states;

    collected_states.clear();
    collected_Y_states_no.clear();
    collected_M_states_no.clear();
    collected_G_states_no.clear();
    collected_killedinIC_states_no.clear();
    collected_antimicrob_effects.clear();
    number_PhagocytosisEvents.clear();

    // epsilon calculation
    relativeSubunitEpsilons.clear();
    bool isEquilibrated = false;

    mek = 0.0;
    mmk = 0.0;
    gek = 0.0;
    ggk = 0.0;
}

void PI_SBM_Solver::simulate(std::minstd_rand& engine) {
    // initialize
    time = 0.0;
    record_time = 0.0;
    curr_sim_step = 0;
    overall_sum_phag = 0.0;
    int counter = 0;
    setModel();
    initializeAntimicrobialEffect();
    initializeResistance();
    int counter_n = 0;

    while (time <= run_time) {
        counter = 0.0;
        counter_n = 0;
        if (record_time == time) { //record_time gives at which time interval we collect the data for the simulation output
            data_collection();
            record_time = record_time + record_data_time_step;
        }
        time = time + time_step;
        curr_sim_step++;

        // variables that must be reset for each time step
        no_phag_l = 0.0;
        no_phag_d = 0.0;
        toInsertG.clear();
        toInsertM.clear();
        toInsertGA.clear();
        toInsertMA.clear();

        // go over all individuals: sum_individuals; don't forget that these might change over time
        sum_Individuals = prelim_YAE + prelim_YKE + start_no_M + start_no_G;
        for (int ind = 0; ind < sum_Individuals; ind++) {
            counter_n++;
            double rand_number = std::uniform_real_distribution<double>{-sum_Individuals,sum_Individuals}(engine);
            if(abs(rand_number)<=prelim_YAE + prelim_YKE){
                updateYstate(abs(rand_number),(rand_number > 0),engine);
                counter++;
            }else if(abs(rand_number)<=prelim_YAE + prelim_YKE + start_no_G){
                updateGstate(abs(rand_number)-prelim_YAE - prelim_YKE,(rand_number > 0),engine);

            }else{
                updateMstate(abs(rand_number)-prelim_YAE - prelim_YKE - start_no_G,(rand_number > 0),engine);
            }
            sum_Individuals = prelim_YAE + prelim_YKE + start_no_M + start_no_G; //
        }

        updateICset(); // IC set and ICA set
        // update IC states!!! --> so that not 2 processes are performed during one time step
        updateICstateVector();

        // calculate relative subunit epsilons and sets
        updateCurrentYStates();  //the Y in variables or function names refers to "Yeast"

        setAntimicrobEffects();
        if(resistance_effect){
            setResistance();
        }
    }
}


void PI_SBM_Solver::updateYstate(double selected_agent, bool reverse, std::minstd_rand& engine) {

    bool reaction = false;

    // select Y out of all Y that can perform a reaction (without YR, YGK, YMK)
    int Ystate = selectYstate(selected_agent); // 1 = YAE, 2=YKE, 3=YGA, 4=YMA

    switch (Ystate) {
    case 1:
        if (reverse) {
            reaction = testForResistance(1,engine);
            if (!reaction) {
                reaction = testForExtracellularKilling(engine);
            } else {
                break;
            }
        } else {
            reaction = testForExtracellularKilling(engine);
            if (!reaction) {
                reaction = testForResistance(1,engine);
            }
        }
        break;
    case 2:
        reaction = testForResistance(2,engine);
        break;
    }
}

void PI_SBM_Solver::updateGstate(double selected_agent,bool reverse, std::minstd_rand& engine) {
    const auto[i, j] = selectIndividualOfICStateG(selected_agent);
    if (i == 0 || (!reverse and !testForIntracellularKillingG(i, j, selected_agent-floor(selected_agent)))) { //most used case
        testForPhagocytosisG(i, j, engine);
    }
    else if (reverse) {
        if (!testForPhagocytosisG(i, j, engine)) {
            testForIntracellularKillingG(i, j, selected_agent-floor(selected_agent));
        }
    }
}
void PI_SBM_Solver::updateMstate(double selected_agent, bool reverse, std::minstd_rand& engine) {
    const auto [i, j] = selectIndividualOfICStateM(selected_agent);
    if (i == 0 || (!reverse and !testForIntracellularKillingM(i, j, selected_agent-floor(selected_agent)))) { //most used case
        testForPhagocytosisM(i, j, engine);
    } else if (reverse) {
        if (!testForPhagocytosisM(i, j, engine)) {
            testForIntracellularKillingM(i, j, selected_agent-floor(selected_agent));
        }
    }
}

std::pair<int,int> PI_SBM_Solver::selectIndividualOfICStateG(double selected_agent) {
    double cell_count = 0;
    if(prelim_G_states[0][0] >= selected_agent){//most used case
        return {0,0};
    }
    for(const auto& temp_xa: occupied_G_states){
        int i = (temp_xa / 10);
        int j = (temp_xa - (i * 10));
        cell_count += prelim_G_states[i][j];
        if(selected_agent <= cell_count){
            return {i,j};
        }
    }
    return {0,0};
}

std::pair<int,int> PI_SBM_Solver::selectIndividualOfICStateM(double selected_agent) {
    double cell_count = 0;
    if(prelim_M_states[0][0] >= selected_agent){//most used case
        return {0,0};
    }
    for(const auto& temp_xa: occupied_M_states){
        int i = (temp_xa / 10);
        int j = (temp_xa - (i * 10));
        cell_count += prelim_M_states[i][j];
        if(selected_agent <= cell_count){
            return {i,j};
        }
    }
    return {0,0};
}

int PI_SBM_Solver::selectYstate(double selected_agent) const {
    int Ystate = 0;
    if (selected_agent <= prelim_YAE) { // YAE
        Ystate = 1;
    } else {
        Ystate = 2;
    }
    return Ystate;
}

bool PI_SBM_Solver::testForExtracellularKilling(std::minstd_rand& engine) {
    bool reaction = false;
    double rand_numb = std::uniform_real_distribution<double>{0.0,1.0}(engine);

    if (rand_numb <= ((collected_antimicrob_effects[curr_sim_step - 1][4]) * time_step)) { //(kappaEK*time_step)){
        prelim_YAE--;
        toAddYKE++;
        reaction = true;
    } else {
        reaction = false;
    }
    return reaction;
}

bool PI_SBM_Solver::testForResistance(int Ystate,std::minstd_rand& engine) {
    bool reaction = false;
    double rand_numb = std::uniform_real_distribution<double>{0.0,1.0}(engine); // generateRandomRealNumber(ub, lb);
    double probability = rho * time_step;

    if (rand_numb <= probability) { //(kappaEK*time_step)){
        reaction = true;
        switch (Ystate) {
        case 1:
            prelim_YAE--;
            toAddYRA++;
            break;
        case 2:
            prelim_YKE--;
            toAddYRK++;
            break;
        }
    } else {
        reaction = false;
    }
    return reaction;
}

bool PI_SBM_Solver::testForIntracellularKillingG(int i, int j, double probability) {
    bool success = false;
    double kappa = 0; //, ub_d;
    double ub_d = no_d - 1;

    int temp_i = 0, temp_j = 0, temp_no = 0, temp_x = 0;


    if (i < 1) {
        success = false;
    } else {
            kappa = kappaG;

        if (probability <= (kappa * time_step)) {

            if (j < ub_d) {

                // perform reaction
                success = true;
                prelim_G_states[i][j] = prelim_G_states[i][j] - 1;

                temp_no = static_cast<int>(prelim_G_states[i][j]); // ICK_current_IC_states.at(i).at(j);

                if (temp_no == 0) {
                    temp_x = getSetValueFromIndizes(i, j);
                    occupied_G_states.erase(occupied_G_states.find(temp_x));
                    occupied_GA_states.erase(occupied_GA_states.find(temp_x));
                }

                temp_i = i - 1;
                temp_j = j + 1;
                prelim_G_states[temp_i][temp_j] = prelim_G_states[temp_i][temp_j] + 1;

                temp_no = static_cast<int>(prelim_G_states[temp_i][temp_j]); // ICK_current_IC_states.at(temp_i).at(temp_j);
                temp_x = getSetValueFromIndizes(temp_i, temp_j);

                if (temp_no == 1) {
                    occupied_G_states.insert(temp_x);
                    if (temp_i > 0) {
                        occupied_GA_states.insert(temp_x);
                    }
                }

                    prelim_YAG = prelim_YAG - 1;
                    ggk = ggk + 1, toAddYKG++;

            }
        } else {
            success = false;
        }
    }
    return success;
}



bool PI_SBM_Solver::testForIntracellularKillingM(int i, int j, double probability) {
    bool success = false;
    double kappa = 0;
    double ub_d = no_d - 1;

    int temp_i = 0, temp_j = 0, temp_no = 0, temp_x = 0;

    if (i < 1) {
        success = false;
    } else {
            kappa = kappaM;

        if (probability <= (kappa * time_step)) { //(kappaEK*time_step)){

            if (j < ub_d) {

                // perform reaction
                success = true;
                prelim_M_states[i][j] = prelim_M_states[i][j] - 1;
                temp_no = static_cast<int>(prelim_M_states[i][j]); // ICK_current_IC_states.at(i).at(j);

                if (temp_no == 0) {
                    temp_x = getSetValueFromIndizes(i, j);
                    occupied_M_states.erase(occupied_M_states.find(temp_x));
                    occupied_MA_states.erase(occupied_MA_states.find(temp_x));
                }

                temp_i = i - 1;
                temp_j = j + 1;
                prelim_M_states[temp_i][temp_j] = prelim_M_states[temp_i][temp_j] + 1;

                temp_no = static_cast<int>(prelim_M_states[temp_i][temp_j]); // ICK_current_IC_states.at(temp_i).at(temp_j);
                temp_x = getSetValueFromIndizes(temp_i, temp_j);
                if (temp_no == 1) {
                    occupied_M_states.insert(temp_x);
                    if (temp_i > 0) {
                        occupied_MA_states.insert(temp_x);
                    }
                }
                prelim_YAM = prelim_YAM - 1;
                mmk = mmk + 1, toAddYKM++;

            }
        } else {
            success = false;
        }
    }
    return success;
}

bool PI_SBM_Solver::testForPhagocytosisM(int i, int j, std::minstd_rand& engine) {
    bool success = false;
    const double rand_numb = std::uniform_real_distribution<double>{0.0,prelim_YAE + prelim_YKE + prelim_YR + start_no_M + start_no_G}(engine);
    int Ystate = 0;
    if (rand_numb <= prelim_YAE) {
        Ystate = 1; // YEA
    }else if (rand_numb <= prelim_YAE + prelim_YKE) {
        Ystate = 2; // YKE
    }
    if (Ystate != 0 && (rand_numb-floor(rand_numb)) <= (phiM * time_step)) {
        double ub_d = no_d - 1;
        double temp_no = 0.0;
        double temp_x = 0.0;
        int temp_i = 0;
        int temp_j = 0;
        int ek = 0;
        if (j < ub_d && i < ub_d) {
            success = true;
            prelim_M_states[i][j] = prelim_M_states[i][j] - 1;
            temp_no = prelim_M_states[i][j];

            if (temp_no == 0) {
                temp_x = getSetValueFromIndizes(i, j);
                occupied_M_states.erase(occupied_M_states.find(temp_x));
                if (i > 0) {
                    occupied_MA_states.erase(occupied_MA_states.find(temp_x));
                }
            }

            switch (Ystate) {
            case 1:
                temp_i = i + 1;
                temp_j = j;
                prelim_YAE--;
                ek = 0;

                break; // phagocytosis of an alive Y
            case 2:
                temp_i = i;
                temp_j = j + 1;
                prelim_YKE--;
                ek = ek + 1;

                break; // phagocytosis of a killed Y
            }
            prelim_M_states[temp_i][temp_j] = prelim_M_states[temp_i][temp_j] + 1;
            temp_no = prelim_M_states[temp_i][temp_j]; // = temp_M.at(rand_l).at(rand_d+1) + 1;
            temp_x = getSetValueFromIndizes(temp_i, temp_j);

            if (temp_no == 1) {
                occupied_M_states.insert(temp_x);
                if (temp_i > 0) {
                    occupied_MA_states.insert(temp_x);
                }
            }

                if (ek == 0) {
                    toAddYAM++;
                } else {
                    mek = mek + ek;
                    toAddYKM++;
                }
        }
    }
    return success;
}

bool PI_SBM_Solver::testForPhagocytosisG(int i, int j, std::minstd_rand& engine) {
    bool success = false;
    double phi;
    if (i == 0 && j == 0) {
        phi = phiG;
    } else {
        phi = phiGstar;
    }
    const double rand_numb = std::uniform_real_distribution<double>{0.0,prelim_YAE + prelim_YKE + prelim_YR + start_no_M + start_no_G}(engine);
    int Ystate = 0;
    if (rand_numb <= prelim_YAE) {
        Ystate = 1; // YEA
    }else if (rand_numb <= prelim_YAE + prelim_YKE) {
        Ystate = 2; // YKE
    }
    if (Ystate != 0 && (rand_numb-floor(rand_numb)) <= (phi * time_step)) {
        double ub_d = no_d - 1;
        double temp_no = 0.0;
        double temp_x = 0.0;
        int temp_i = 0;
        int temp_j = 0;
        int ek = 0;
        if (j < ub_d && i < ub_d) {
            success = true;
            prelim_G_states[i][j] = prelim_G_states[i][j] - 1;
            temp_no = prelim_G_states[i][j];

            if (temp_no == 0) {
                temp_x = getSetValueFromIndizes(i, j);
                occupied_G_states.erase(occupied_G_states.find(temp_x));
                if (i > 0) {
                    occupied_GA_states.erase(occupied_GA_states.find(temp_x));
                }
            }

            switch (Ystate) {
            case 1:
                temp_i = i + 1;
                temp_j = j;
                prelim_YAE--;
                ek = 0;
                if (i == 0 && j == 0) {
                    no_phag_l++;
                };
                break; // phagocytosis of an alive Y
            case 2:
                temp_i = i;
                temp_j = j + 1;
                prelim_YKE--;
                ek = ek + 1;
                if (i == 0 && j == 0) {
                    no_phag_d++;
                };
                break; // phagocytosis of a killed Y
            }
            prelim_G_states[temp_i][temp_j] = prelim_G_states[temp_i][temp_j] + 1;
            temp_no = prelim_G_states[temp_i][temp_j];
            temp_x = getSetValueFromIndizes(temp_i, temp_j);

            if (temp_no == 1) {
                occupied_G_states.insert(temp_x);
                if (temp_i > 0) {
                    occupied_GA_states.insert(temp_x);
                }
            }

            if (ek == 0) {
                toAddYAG++;
            } else {
                gek = gek + ek;
                toAddYKG++;
            }

        }
    }
    return success;
}

void PI_SBM_Solver::initializeResistance(){
    collected_resistance_effects.clear();
    std::vector<double> temp_resistance_eff(5,0.0);
    collected_resistance_effects.push_back(temp_resistance_eff);

}

void PI_SBM_Solver::initializeAntimicrobialEffect() {

    collected_antimicrob_effects.clear();

    collected_no_KE.clear();
    number_PhagocytosisEvents.push_back(0.0);
    std::vector<double> temp_antimi_eff(5, 0.0);
    collected_antimicrob_effects.push_back(temp_antimi_eff);
    std::vector<double> temp_no_KE(4, 0.0);
    collected_no_KE.push_back(temp_no_KE);
}
void PI_SBM_Solver::setAntimicrobEffects() {

    double sum = 0.0, G10_eff = 0.0, G01_eff = 0.0, kek = 0.0;
    sum = 0.0;
    kek = 0.0;

    G10_eff = no_phag_l; /// start_no_G;
    G01_eff = no_phag_d; /// start_no_G;
    sum = G10_eff + G01_eff;
    number_PhagocytosisEvents.push_back(sum);

    // calc new effect:
    double tt = 0.0;
    int counter = 0;
    while (tt <= time) {

        overall_sum_phag = overall_sum_phag + (number_PhagocytosisEvents[counter] * exp(-gamma * (time - tt)));
        counter++;
        tt = tt + time_step;
    }
    overall_sum_phag = overall_sum_phag * kappaEG / start_no_G;
    std::vector<double> temp_antimi_eff;// {time, sum, gamma, kappaEG, overall_sum_phag};
    temp_antimi_eff.push_back(time);
    temp_antimi_eff.push_back(sum); // phag events
    temp_antimi_eff.push_back(gamma);
    temp_antimi_eff.push_back(kappaEG);
    temp_antimi_eff.push_back(overall_sum_phag);
    collected_antimicrob_effects.push_back(temp_antimi_eff);
}
/**
 * calculate new resistance rate dependend on phagocytosis events and time
 */
void PI_SBM_Solver::setResistance()
{
    //calc new effect:
    double tt = 0.0;
    int counter = 0;
    while( tt <= time )
    {
        overall_sum_resistance = overall_sum_resistance + ( number_PhagocytosisEvents[counter]*exp(-gammaR*time_step*( time - tt ) ) );
        counter++;
        tt = tt+time_step;
    }
    overall_sum_resistance = overall_sum_resistance*rho / start_no_G;// start_no_G; //
    std::vector<double> temp_resistance;// {time, sum, gammaR, rhobar, overall_sum_resistance};
    temp_resistance.push_back(time);
    temp_resistance.push_back(no_phag_l + no_phag_d); //phag events
    temp_resistance.push_back(gammaR);
    temp_resistance.push_back(rho);
    temp_resistance.push_back(overall_sum_resistance);
    collected_resistance_effects.push_back(temp_resistance);
}

void PI_SBM_Solver::updateICstateVector() {
    temp_G_states = prelim_G_states;
    temp_M_states = prelim_M_states;
}

void PI_SBM_Solver::updateICset() {
    std::set<int>::iterator it_G2;
    for (int i = 0; i < toInsertG.size(); i++) {
        it_G2 = find(occupied_G_states.begin(), occupied_G_states.end(), toInsertG[i]);
        if (it_G2 != occupied_G_states.end()) {
        } else {
            occupied_G_states.insert(toInsertG[i]);
        }
    }
    std::set<int>::iterator it_M2;
    for (int i = 0; i < toInsertM.size(); i++) {
        it_M2 = find(occupied_M_states.begin(), occupied_M_states.end(), toInsertM[i]);
        if (it_M2 != occupied_M_states.end()) {
        } else {
            occupied_M_states.insert(toInsertM[i]);
        }
    }

    std::set<int>::iterator it_GA2;
    for (int i = 0; i < toInsertGA.size(); i++) {
        it_GA2 = find(occupied_GA_states.begin(), occupied_GA_states.end(), toInsertGA[i]);
        if (it_GA2 != occupied_GA_states.end()) {
        } else {
            occupied_GA_states.insert(toInsertGA[i]);
        }
    }

    std::set<int>::iterator it_MA2;
    for (int i = 0; i < toInsertMA.size(); i++) {
        it_MA2 = find(occupied_MA_states.begin(), occupied_MA_states.end(), toInsertMA[i]);
        if (it_MA2 != occupied_MA_states.end()) {
        } else {
            occupied_MA_states.insert(toInsertMA[i]);
        }
    }
}

void PI_SBM_Solver::updateCurrentYStates() {

    // create a temporary vector of epsilon errors
    std::vector<double> temp(9);
    temp[0] = time;
    temp[1] = abs((prelim_YAE + toAddYAE - temp_YAE) / temp_YAE);
    temp[2] = abs((prelim_YKE + toAddYKE - temp_YKE) / temp_YKE);
    temp[3] = abs((prelim_YRA + toAddYRA - temp_YRA) / temp_YRA);
    temp[4] = abs((prelim_YRK + toAddYRK - temp_YRK) / temp_YRK);
    temp[5] = abs((prelim_YAM + toAddYAM - temp_YAM) / temp_YAM);
    temp[6] = abs((prelim_YAG + toAddYAG - temp_YAG) / temp_YAG);
    temp[7] = abs((prelim_YKM + toAddYKM - temp_YKM) / temp_YKM);
    temp[8] = abs((prelim_YKG + toAddYKG - temp_YKG) / temp_YKG);

    // append vector to epsilon matrix
    relativeSubunitEpsilons.push_back(temp);

    temp_YAE = prelim_YAE + toAddYAE;
    temp_YKE = prelim_YKE + toAddYKE;

    temp_YRA = prelim_YRA + toAddYRA;
    temp_YRK = prelim_YRK + toAddYRK;
    temp_YR = temp_YRA + temp_YRK;
    temp_YAM = prelim_YAM + toAddYAM;
    temp_YAG = prelim_YAG + toAddYAG;
    temp_YKM = prelim_YKM + toAddYKM;
    temp_YKG = prelim_YKG + toAddYKG;

    // reset for next timestep:
    prelim_YAE = temp_YAE;
    prelim_YKE = temp_YKE;
    prelim_YR = temp_YR;
    prelim_YRA = temp_YRA;
    prelim_YRK = temp_YRK;

    prelim_YAM = temp_YAM;
    prelim_YAG = temp_YAG;
    prelim_YKM = temp_YKM;
    prelim_YKG = temp_YKG;

    toAddYAE = 0.0;
    toAddYKE = 0.0;
    toAddYAG = 0.0;
    toAddYAM = 0.0;
    toAddYRA = 0.0;
    toAddYRK = 0.0;
    toAddYKG = 0.0;
    toAddYKM = 0.0;
}


void PI_SBM_Solver::data_collection() {

    std::vector<double> temp_states;
    std::vector<double> temp_states_Y;
    std::vector<double> temp_states_M_j;
    std::vector<std::vector<double>> temp_states_M;
    std::vector<double> temp_states_G_j;
    std::vector<std::vector<double>> temp_states_G;

    temp_states.push_back(time);
    temp_states.push_back((double)(temp_YAE / temp_Y0));
    temp_states.push_back((double)(temp_YR / temp_Y0));
    temp_states.push_back((double)(temp_YKE / temp_Y0));

    temp_states_Y.push_back(time);
    temp_states_Y.push_back((double)(temp_YAE));
    temp_states_Y.push_back((double)(temp_YR));
    temp_states_Y.push_back((double)(temp_YKE));

    temp_states_M.clear();
    for (int i = 0; i < temp_M_states.size(); i++) {

        temp_states_M_j.clear();
        for (int j = 0; j < temp_M_states[i].size(); j++) {
            temp_states.push_back(temp_M_states[i][j] / start_no_M);
            temp_states_M_j.push_back(temp_M_states[i][j]);
        }
        temp_states_M.push_back(temp_states_M_j);
    }
    temp_states_G.clear();
    for (int i = 0; i < temp_G_states.size(); i++) {
        temp_states_G_j.clear();
        for (int j = 0; j < temp_G_states[i].size(); j++) {
            temp_states.push_back(temp_G_states[i][j] / start_no_G);
            temp_states_G_j.push_back(temp_G_states[i][j]);
        }
        temp_states_G.push_back(temp_states_G_j);
    }

    collected_states.push_back(temp_states);
    collected_Y_states_no.push_back(temp_states_Y);
    collected_M_states_no.push_back(temp_states_M);
    collected_G_states_no.push_back(temp_states_G);

    std::vector<double> temp_collected_killedinIC_states_no; // {time, (mmk/temp_Y0), (ggk/temp_Y0), (mek/temp_Y0), (gek/temp_Y0)};
    temp_collected_killedinIC_states_no.push_back(time);
    temp_collected_killedinIC_states_no.push_back(mmk / temp_Y0);
    temp_collected_killedinIC_states_no.push_back(ggk / temp_Y0);
    temp_collected_killedinIC_states_no.push_back(mek / temp_Y0);
    temp_collected_killedinIC_states_no.push_back(gek / temp_Y0);
    collected_killedinIC_states_no.push_back(temp_collected_killedinIC_states_no);

    std::vector<double> temp_collected_KE; // {time, (temp_YRA/temp_Y0), (temp_YRK/temp_Y0), (temp_YR/temp_Y0)};
    temp_collected_KE.push_back(time);
    temp_collected_KE.push_back(temp_YRA / temp_Y0);
    temp_collected_KE.push_back(temp_YRK / temp_Y0);
    temp_collected_KE.push_back(temp_YR / temp_Y0);
    collected_no_KE.push_back(temp_collected_KE);
}

int PI_SBM_Solver::getSetValueFromIndizes(int i, int j){
    int x = (i * 10) + j;
    return x;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_MY() {

    double time, sumMY, temp_MYY_val, temp_MYM_val, temp_MYKM, temp_M, temp_MYKE;

    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_my{};

    for (int z = 0; z < collected_M_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        sumMY = 0.0;
        temp_MYY_val = 0.0;

        for (int i = 0; i < collected_M_states_no.at(z).size(); i++) {           // alive
            for (int j = 0; j < collected_M_states_no.at(z).at(i).size(); j++) { // dead
                sumMY = sumMY + (collected_M_states_no.at(z).at(i).at(j) * (i + j));
            }
        }

        temp_MYY_val = sumMY / temp_Y0;
        temp_v.push_back(temp_MYY_val);
        values_sim_data_my.push_back(temp_v);

        temp_v.clear();

    }
    return values_sim_data_my;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_GY() {

    double time, sumGY, temp_GYY_val, temp_GYG_val, temp_G_val, temp_GC5A, temp_GS, temp_GSYKG, temp_GSYKE;

    std::vector<double> temp_v;

    std::vector<std::vector<double>> values_sim_data_gy{};
    for (int z = 0; z < collected_G_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        sumGY = 0.0;
        temp_GYY_val = 0.0;

        for (int i = 0; i < collected_G_states_no.at(z).size(); i++) {           // alive
            for (int j = 0; j < collected_G_states_no.at(z).at(i).size(); j++) { // dead
                sumGY = sumGY + (collected_G_states_no.at(z).at(i).at(j) * (i + j));
            }
        }

        temp_GYY_val = sumGY / temp_Y0;
        temp_v.push_back(temp_GYY_val);
        values_sim_data_gy.push_back(temp_v);

        temp_v.clear();
    }
    return values_sim_data_gy;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YGA() {

    double time, sumYGA, temp_YGA_val;

    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_yga{};

    for (int z = 0; z < collected_G_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);

        sumYGA = 0.0;
        temp_YGA_val = 0.0;
        auto alive = 0;
        for (int i = 0; i < collected_G_states_no.at(z).size(); i++) {
            for (int j = 0; j < collected_G_states_no.at(z).at(i).size(); j++) {
                sumYGA = sumYGA + (collected_G_states_no.at(z).at(i).at(j) * i);  // alive
            }
        }

        temp_YGA_val = sumYGA / temp_Y0;
        temp_v.push_back(temp_YGA_val);
        values_sim_data_yga.push_back(temp_v);

        temp_v.clear();
    }
    return values_sim_data_yga;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YMA() {

    double time, sumYMA, temp_YMA_val;

    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_yma{};

    for (int z = 0; z < collected_M_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);

        sumYMA = 0.0;
        temp_YMA_val = 0.0;
        auto alive = 0;
        for (int i = 0; i < collected_M_states_no.at(z).size(); i++) {
            for (int j = 0; j < collected_M_states_no.at(z).at(i).size(); j++) {
                sumYMA = sumYMA + (collected_M_states_no.at(z).at(i).at(j) * i);  // alive
            }
        }

        temp_YMA_val = sumYMA / temp_Y0;
        temp_v.push_back(temp_YMA_val);
        values_sim_data_yma.push_back(temp_v);

        temp_v.clear();
    }
    return values_sim_data_yma;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YF() {
    double time, temp_YH_val, temp_Y_val, temp_YKE_val, temp_YF_val;
    std::vector<double> temp_v;

    std::vector<std::vector<double>> values_sim_data_yf{};

    for (int i = 0; i < collected_Y_states_no.size(); i++) {
        time = collected_Y_states_no.at(i).at(0);
        temp_v.push_back(time);

        temp_Y_val = collected_Y_states_no.at(i).at(1) / temp_Y0;
        temp_YH_val = collected_Y_states_no.at(i).at(2) / temp_Y0;
        temp_YKE_val = collected_Y_states_no.at(i).at(3) / temp_Y0;
        temp_YF_val = temp_YH_val + temp_Y_val + temp_YKE_val;

        temp_v.push_back(temp_YF_val);
        temp_v.push_back(temp_Y_val);
        temp_v.push_back(temp_YH_val);
        temp_v.push_back(temp_YKE_val);

        values_sim_data_yf.push_back(temp_v);

        temp_v.clear();

    }
    return values_sim_data_yf;
}


std::vector<std::vector<double>> PI_SBM_Solver::calculate_YL() {
    double time, temp_MY_val, temp_GY_val, temp_Y_val, temp_YH_val, temp_YL_val;
    std::vector<double> temp_v;

    std::vector<std::vector<double>>  values_sim_data_yl{};

    for (int z = 0; z < collected_Y_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);

        temp_Y_val = collected_Y_states_no.at(z).at(1) / temp_Y0;
        temp_YH_val = collected_no_KE.at(z + 1).at(1); // only alive resistant Yra ----- sim_Y_states.at(z).at(2)/Y_0;at z-1, as KE in SBM model should be initialized


        temp_GY_val = 0.0;
        for (int i = 1; i < collected_G_states_no.at(z).size(); i++) {           // alive
            for (int j = 0; j < collected_G_states_no.at(z).at(i).size(); j++) { // dead
                temp_GY_val = temp_GY_val + (collected_G_states_no.at(z).at(i).at(j) * i);
            }
        }
        temp_MY_val = 0.0;
        for (int i = 1; i < collected_M_states_no.at(z).size(); i++) {           // alive
            for (int j = 0; j < collected_M_states_no.at(z).at(i).size(); j++) { // dead
                temp_MY_val = temp_MY_val + (collected_M_states_no.at(z).at(i).at(j) * i);
            }
        }
        temp_GY_val = temp_GY_val / temp_Y0;
        temp_MY_val = temp_MY_val / temp_Y0;
        temp_YL_val = temp_Y_val + temp_YH_val + temp_GY_val + temp_MY_val;
        temp_v.push_back(temp_YL_val);
        temp_v.push_back(temp_Y_val);
        temp_v.push_back(temp_YH_val);
        temp_v.push_back(temp_MY_val);
        temp_v.push_back(temp_GY_val);

        values_sim_data_yl.push_back(temp_v);

        temp_v.clear();
    }
    return values_sim_data_yl;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YK() {
    double time, temp_YKE_val, temp_MY_val, temp_GY_val, temp_YK_val;

    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_YK{};

    for (int z = 0; z < collected_Y_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_YKE_val = (collected_Y_states_no.at(z).at(3) / temp_Y0) + collected_no_KE.at(z + 1).at(2); // add Yrk at z-1
        temp_GY_val = 0.0;

        for (int i = 0; i < collected_G_states_no.at(z).size(); i++) {           // alive
            for (int j = 1; j < collected_G_states_no.at(z).at(i).size(); j++) { // dead
                temp_GY_val = temp_GY_val + (collected_G_states_no.at(z).at(i).at(j) * j);
            }
        }
        temp_MY_val = 0.0;
        for (int i = 0; i < collected_M_states_no.at(z).size(); i++) {           // alive
            for (int j = 1; j < collected_M_states_no.at(z).at(i).size(); j++) { // dead
                temp_MY_val = temp_MY_val + (collected_M_states_no.at(z).at(i).at(j) * j);
            }
        }
        temp_MY_val = temp_MY_val / temp_Y0;
        temp_GY_val = temp_GY_val / temp_Y0;
        temp_YK_val = temp_YKE_val + temp_MY_val + temp_GY_val;

        temp_v.push_back(temp_YK_val);
        temp_v.push_back(temp_YKE_val);
        temp_v.push_back(temp_GY_val);
        temp_v.push_back(temp_MY_val);

        values_sim_data_YK.push_back(temp_v);

        temp_v.clear();
    }
    return values_sim_data_YK;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YKG(){
    double time, temp_GY_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_YKG{};
    for (int z = 0; z < collected_Y_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_GY_val = 0.0;
        for (int i = 0; i < collected_G_states_no.at(z).size(); ++i) {           // lebende
            for (int j = 1; j < collected_G_states_no.at(z).at(i).size(); ++j) { // tote
                temp_GY_val += (collected_G_states_no.at(z).at(i).at(j) * j);
            }
        }
        temp_GY_val = temp_GY_val / temp_Y0;
        temp_v.push_back(temp_GY_val);
        values_sim_data_YKG.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_YKG;
}
std::vector<std::vector<double>> PI_SBM_Solver::calculate_YKM(){
    double time, temp_MY_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_YKM{};
    for (int z = 0; z < collected_Y_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_MY_val = 0.0;
        for (int i = 0; i < collected_M_states_no.at(z).size(); i++) {           // lebende
            for (int j = 1; j < collected_M_states_no.at(z).at(i).size(); j++) { // tote
                temp_MY_val += (collected_M_states_no.at(z).at(i).at(j) * j);
            }
        }
        temp_MY_val = temp_MY_val / temp_Y0;
        temp_v.push_back(temp_MY_val);
        values_sim_data_YKM.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_YKM;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YRK(){
    // immune-evasive killed
    double time, temp_YRK_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_YRK{};
    for (int z = 0; z < collected_Y_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_YRK_val = collected_no_KE.at(z+1).at(2);
        temp_v.push_back(temp_YRK_val);
        values_sim_data_YRK.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_YRK;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YRA(){
    // immune-evasive alive
    double time, temp_YRA_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_YRA{};
    for (int z = 0; z < collected_Y_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_YRA_val = collected_no_KE.at(z+1).at(1);
        temp_v.push_back(temp_YRA_val);
        values_sim_data_YRA.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_YRA;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YKE(){
    // killed extracellular (no killed immune-evasive)
    double time, temp_YKE_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_YKE{};
    for (int z = 0; z < collected_Y_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        //temp_YKE_val = temp_YKE/temp_Y0 + temp_YRK/temp_Y0
        temp_YKE_val = (collected_Y_states_no.at(z).at(3) / temp_Y0);
        temp_v.push_back(temp_YKE_val);
        values_sim_data_YKE.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_YKE;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YAE(){
    // alive extracellular (no killed immune-evasive)
    double time, temp_YAE_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_YAE{};
    for (int z = 0; z < collected_Y_states_no.size(); z++) {
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_YAE_val = (collected_Y_states_no.at(z).at(1) / temp_Y0);
        temp_v.push_back(temp_YAE_val);
        values_sim_data_YAE.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_YAE;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_MMK(){
    // phagocytosed by monocytes and killed intracellularly
    double time, temp_MMK_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_MMK{};
    for (int z = 0; z < collected_killedinIC_states_no.size(); z++){
        time = collected_killedinIC_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_MMK_val = (collected_killedinIC_states_no.at(z).at(1));
        temp_v.push_back(temp_MMK_val);
        values_sim_data_MMK.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_MMK;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_GGK(){
    // phagocytosed by monocytes and killed intracellularly
    double time, temp_GGK_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_GGK{};
    for (int z = 0; z < collected_killedinIC_states_no.size(); z++){
        time = collected_killedinIC_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_GGK_val = (collected_killedinIC_states_no.at(z).at(2));
        temp_v.push_back(temp_GGK_val);
        values_sim_data_GGK.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_GGK;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_MEK(){
    // phagocytosed by monocytes and killed intracellularly
    double time, temp_MEK_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_MEK{};
    for (int z = 0; z < collected_killedinIC_states_no.size(); z++){
        time = collected_killedinIC_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_MEK_val = (collected_killedinIC_states_no.at(z).at(3));
        temp_v.push_back(temp_MEK_val);
        values_sim_data_MEK.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_MEK;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_GEK(){
    // phagocytosed by monocytes and killed intracellularly
    double time, temp_GEK_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_GEK{};
    for (int z = 0; z < collected_killedinIC_states_no.size(); z++){
        time = collected_killedinIC_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_GEK_val = (collected_killedinIC_states_no.at(z).at(4));
        temp_v.push_back(temp_GEK_val);
        values_sim_data_GEK.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_GEK;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_PE(){ //# phagocytosis events
    // phagocytosed by monocytes and killed intracellularly
    double time, temp_PE_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_PE{};
    for (int z = 0; z < collected_killedinIC_states_no.size(); z++){
        time = collected_killedinIC_states_no.at(z).at(0);
        temp_v.push_back(time);
        temp_PE_val = number_PhagocytosisEvents.at(z) / start_no_G;
        temp_v.push_back(temp_PE_val);
        values_sim_data_PE.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_PE;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_kEK(){
    //calculate the time-varying extracellular killing rate
    double time, temp_kEK_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_kEK{};
    for (int z = 0; z < collected_antimicrob_effects.size(); z++){
        time = collected_antimicrob_effects.at(z).at(0);
        temp_v.push_back(time);
        temp_kEK_val = collected_antimicrob_effects.at(z).at(4);
        temp_v.push_back(temp_kEK_val);
        values_sim_data_kEK.push_back(temp_v);
        temp_v.clear();
    }
    return values_sim_data_kEK;
}

std::vector<std::vector<double>> PI_SBM_Solver::calculate_YR(){
    // calculate the number of resistant pathogens
    double time, temp_YR_val;
    std::vector<double> temp_v;
    std::vector<std::vector<double>> values_sim_data_YR{};
    for (int z = 0; z < collected_Y_states_no.size(); z++){
        time = collected_Y_states_no.at(z).at(0);
        temp_v.push_back(time);
            temp_YR_val = collected_no_KE.at(z+1).at(2) + collected_no_KE.at(z+1).at(1);
            temp_v.push_back(temp_YR_val);
            values_sim_data_YR.push_back(temp_v);
            temp_v.clear();
    }
    return values_sim_data_YR;
}
