/*
 * File:   main.cpp
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

#include "macros.h"
#include <boost/filesystem.hpp>
#include "utils/misc.h"
#include "utils/io_utils.h"
#include "model/PI_SBM_Solver.h"
#include <fstream>
#include <chrono>
#include <random>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)



int main(int argc, char** argv) {
    //load parameters from config file (json)
    boost::filesystem::path config_json("../../config.json");
    if (argc <= 1) {
        ERROR_STDERR("usage: " << argv[0] << " <absolute_path_to_json>");
        return 1;
    }
    if (!(std::istringstream(argv[1]) >> config_json)) {
        ERROR_STDERR("usage: " << argv[0] << " <absolute_path_to_json>");
        return 2;
    }


    auto parameters = sbm::util::get_main_config_parameters(config_json.string());

    std::minstd_rand engine{parameters.general_parameters.standard_seed};
/*
 * Initial Values:
    parameters.model_parameters.cells["ALIVEEXTERNCANDIDA"];  // Alive Extracelluar Candida
    parameters.model_parameters.cells["RESISTENTCANDIDA"];    // Total Resistent Extracelluar Candida
    parameters.model_parameters.cells["ALIVEIECANDIDA"];      // Alive Resistent Extracelluar Candida
    parameters.model_parameters.cells["KILLEDIECANDIDA"];     // Killed Resistent Extracelluar Candida
    parameters.model_parameters.cells["KILLEDEXTERNCANDIDA"]; // Killed Extracelluar Candida
    parameters.model_parameters.cells["MONOCYTES"];           // Number of Monocytes
    parameters.model_parameters.cells["GRANULOCYTES"];      // Number of Granulocytes

    Initial Rates, Range [0,1]:
    parameters.model_parameters.rates["PhiG"];     // Phagocytosis Granulocytes (First)
    parameters.model_parameters.rates["PhiGStar"]; // Phagocytosis Granulocytes (Consecutive)
    parameters.model_parameters.rates["PhiM"];     // Phagocytosis Monocytes
    parameters.model_parameters.rates["KappaM"];   // Intracelluar Killing Monocytes
    parameters.model_parameters.rates["KappaG"];   // Intracelluar Killing Granulocytes
    parameters.model_parameters.rates["KappaExt"];   // ExtracelluarKilling Candida ALbicans
    parameters.model_parameters.rates["Rho"];      // ResistentRate Candida Albicans
    parameters.model_parameters.rates["Gamma"];    // ExtracelluarKilling Candida Albicans(Antimicrobial effect)
    parameters.model_parameters.rates["GammaR"];    // Resistent depending on Neutrophil-phagocyted Candida Albicans(Resistent effect)

    Set parameters like this:
    parameters.model_parameters.rates["Gamma"] = 0.5;
    */
    const auto parameters_simulation = sbm::util::generate_parameter_vector(parameters.model_parameters);

    using timer = std::chrono::steady_clock;
    const auto begin = timer::now();
    const auto solver = std::make_unique<PI_SBM_Solver>(parameters.simulation_parameters,parameters_simulation, parameters.model_parameters.resistance_effect);

    solver->simulate(engine);

    std::ifstream infile;
    infile.open(parameters.general_parameters.input_folder+"/whole_expData_adapted.dat", std::ifstream::in);
    double ca_f_value{};
    double ca_m_ca_value{};
    double ca_g_ca_value{};
    double ca_k_value{};
    double ca_ae_value{};
    double _{};
    double time{};
    double squared_error = 0.0;
    auto values_sim_data_yk = solver->calculate_YK();
    auto values_sim_data_yf = solver->calculate_YF();
    auto values_sim_data_my = solver->calculate_MY();
    auto values_sim_data_gy = solver->calculate_GY();
    auto values_sim_data_yl = solver->calculate_YL();
    auto values_sim_data_ykm = solver->calculate_YKM();
    auto values_sim_data_ykg = solver->calculate_YKG();
    auto values_sim_data_yke = solver->calculate_YKE();
    auto values_sim_data_yae = solver->calculate_YAE();
    auto values_sim_data_yrk = solver->calculate_YRK();
    auto values_sim_data_yra = solver->calculate_YRA();
    auto values_sim_data_mmk = solver->calculate_MMK();
    auto values_sim_data_ggk = solver->calculate_GGK();
    auto values_sim_data_mek = solver->calculate_MEK();
    auto values_sim_data_gek = solver->calculate_GEK();
    auto values_sim_data_yga = solver->calculate_YGA();
    auto values_sim_data_yma = solver->calculate_YMA();
    auto values_sim_data_pe = solver->calculate_PE(); // # phagocytosis events
    auto values_sim_data_kek = solver->calculate_kEK();
    auto values_sim_data_yr = solver->calculate_YR();

    if (infile.good() && infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            iss >> time >> ca_f_value >> _ >> ca_m_ca_value >> _ >> ca_g_ca_value >> _ >> _ >> _ >> _ >> _ >> ca_k_value >> _;
            const auto time_int = static_cast<uint64_t>(round(time));
            squared_error += (values_sim_data_yk.at(time_int).at(1) - ca_k_value) * (values_sim_data_yk.at(time_int).at(1) - ca_k_value);
            squared_error += (values_sim_data_yf.at(time_int).at(1) - ca_f_value) * (values_sim_data_yf.at(time_int).at(1) - ca_f_value);
            squared_error += (values_sim_data_my.at(time_int).at(1) - ca_m_ca_value) * (values_sim_data_my.at(time_int).at(1) - ca_m_ca_value);
            squared_error += (values_sim_data_gy.at(time_int).at(1) - ca_g_ca_value) * (values_sim_data_gy.at(time_int).at(1) - ca_g_ca_value);
        }
    }
    infile.close();

    infile.open(parameters.general_parameters.input_folder+"/CaMG_expData.dat", std::ifstream::in);
    if (infile.good() && infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            iss >> time >> ca_ae_value >> _;
            const auto time_int = static_cast<uint64_t>(round(time));
            squared_error += (values_sim_data_yl.at(time_int).at(1) - ca_ae_value) * (values_sim_data_yl.at(time_int).at(1) - ca_ae_value);
        }
    }
    infile.close();
    const auto end = timer::now();
    std::cout << "Time needed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << '\n';

    std::cout << squared_error << '\n';
    return 0;
}
