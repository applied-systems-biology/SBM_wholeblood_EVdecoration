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
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;



py::handle  run_sbm(int number_simulations,const double time_step, const int max_steps, const py::dict& cell_count,const py::dict& rates, const bool resistance_effect) {
    //load parameters from config file (json)
    sbm::util::Parameters parameters{};

    parameters.simulation_parameters.timestep = time_step;

    parameters.simulation_parameters.timesteps = max_steps;
    parameters.model_parameters.resistance_effect = resistance_effect;
    parameters.model_parameters.cells["ALIVEEXTERNCANDIDA"] = 1e+05;  // Alive Extracelluar Candida
    parameters.model_parameters.cells["RESISTENTCANDIDA"] = 0;    // Total Resistent Extracelluar Candida
    parameters.model_parameters.cells["ALIVEIECANDIDA"] = 0;      // Alive Resistent Extracelluar Candida
    parameters.model_parameters.cells["KILLEDIECANDIDA"] = 0;     // Killed Resistent Extracelluar Candida
    parameters.model_parameters.cells["KILLEDEXTERNCANDIDA"] = 0; // Killed Extracelluar Candida
    parameters.model_parameters.cells["MONOCYTES"] = 50000;           // Number of Monocytes
    parameters.model_parameters.cells["GRANULOCYTES"] = 500000;      // Number of Granulocytes

    parameters.model_parameters.rates["PhiG"] = 0.0321981;     // Phagocytosis Granulocytes (First)
    parameters.model_parameters.rates["PhiGStar"] = 0.0303716; // Phagocytosis Granulocytes (Consecutive)
    parameters.model_parameters.rates["PhiM"] = .0137036;     // Phagocytosis Monocytes
    parameters.model_parameters.rates["KappaM"] = 0.501573;   // Intracelluar Killing Monocytes
    parameters.model_parameters.rates["KappaG"] = 0.0467262;   // Intracelluar Killing Granulocytes
    parameters.model_parameters.rates["Rho"] = 0.00459136;   // ExtracelluarKilling Candida ALbicans
    parameters.model_parameters.rates["Gamma"] = 0.0093675;      // ResistentRate Candida Albicans
    parameters.model_parameters.rates["KappaExt"] = 0.176844;    // ExtracelluarKilling Candida Albicans(Antimicrobial effect)
    parameters.model_parameters.rates["GammaR"] = .0;    // Resistent depending on Neutrophil-phagocyted Candida Albicans(Resistent effect)

    for (const auto& rate : rates) {
        parameters.model_parameters.rates[std::string(py::str(rate.first))] = std::stod(std::string((py::str(rate.second))));
    }

    for (const auto& cell : cell_count) {
        parameters.model_parameters.cells[std::string(py::str(cell.first))] = std::stod(std::string((py::str(cell.second))));
    }
    auto const seed = std::random_device()();
    std::minstd_rand engine{seed};


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
    std::vector<double> output{};
    int num_outputs = 21;

    output.resize(max_steps*num_outputs*number_simulations);
    for(int num=0; num < number_simulations; num++){
        const auto parameters_simulation = sbm::util::generate_parameter_vector(parameters.model_parameters);

        const auto solver = std::make_unique<PI_SBM_Solver>(parameters.simulation_parameters,parameters_simulation, parameters.model_parameters.resistance_effect);
        solver->simulate(engine);
        auto values_sim_data_yk = solver->calculate_YK();
        auto values_sim_data_yf = solver->calculate_YF();
        auto values_sim_data_my = solver->calculate_MY();
        auto values_sim_data_gy = solver->calculate_GY();
        auto values_sim_data_yl = solver->calculate_YL();
        auto values_sim_data_ykm = solver->calculate_YKM();
        auto values_sim_data_ykg = solver->calculate_YKG();
        auto values_sim_data_yke = solver->calculate_YKE();
        auto values_sim_data_yae = solver->calculate_YAE();
        auto values_sim_data_yrk = solver ->calculate_YRK();
        auto values_sim_data_yra = solver->calculate_YRA();
        auto values_sim_data_mmk = solver->calculate_MMK();
        auto values_sim_data_ggk = solver->calculate_GGK();
        auto values_sim_data_mek = solver->calculate_MEK();
        auto values_sim_data_gek = solver->calculate_GEK();
        auto values_sim_data_yga = solver->calculate_YGA();
        auto values_sim_data_yma = solver->calculate_YMA();
        auto values_sim_data_pe = solver->calculate_PE();
        auto values_sim_data_kek = solver->calculate_kEK();
        auto values_sim_data_yr = solver->calculate_YR();

        for(int time_int = 0; time_int < max_steps; ++time_int){
            output[num*num_outputs*max_steps+0+time_int*num_outputs] = values_sim_data_yl.at(time_int).at(1);
            output[num*num_outputs*max_steps+1+time_int*num_outputs] = values_sim_data_yk.at(time_int).at(1);
            output[num*num_outputs*max_steps+2+time_int*num_outputs] = values_sim_data_yf.at(time_int).at(1);
            output[num*num_outputs*max_steps+3+time_int*num_outputs] = values_sim_data_my.at(time_int).at(1);
            output[num*num_outputs*max_steps+4+time_int*num_outputs] = values_sim_data_gy.at(time_int).at(1);
            output[num*num_outputs*max_steps+5+time_int*num_outputs] = values_sim_data_ykm.at(time_int).at(1);
            output[num*num_outputs*max_steps+6+time_int*num_outputs] = values_sim_data_ykg.at(time_int).at(1);
            output[num*num_outputs*max_steps+7+time_int*num_outputs] = values_sim_data_yke.at(time_int).at(1);
            output[num*num_outputs*max_steps+8+time_int*num_outputs] = values_sim_data_yae.at(time_int).at(1);
            output[num*num_outputs*max_steps+9+time_int*num_outputs] = values_sim_data_yrk.at(time_int).at(1);
            output[num*num_outputs*max_steps+10+time_int*num_outputs] = values_sim_data_yra.at(time_int).at(1);
            output[num*num_outputs*max_steps+11+time_int*num_outputs] = values_sim_data_mmk.at(time_int).at(1);
            output[num*num_outputs*max_steps+12+time_int*num_outputs] = values_sim_data_ggk.at(time_int).at(1);
            output[num*num_outputs*max_steps+13+time_int*num_outputs] = values_sim_data_mek.at(time_int).at(1);
            output[num*num_outputs*max_steps+14+time_int*num_outputs] = values_sim_data_gek.at(time_int).at(1);
            output[num*num_outputs*max_steps+15+time_int*num_outputs] = values_sim_data_yga.at(time_int).at(1);
            output[num*num_outputs*max_steps+16+time_int*num_outputs] = values_sim_data_yma.at(time_int).at(1);
            output[num*num_outputs*max_steps+17+time_int*num_outputs] = values_sim_data_yr.at(time_int).at(1);
            output[num*num_outputs*max_steps+18+time_int*num_outputs] = values_sim_data_pe.at(time_int).at(1);
            output[num*num_outputs*max_steps+19+time_int*num_outputs] = values_sim_data_kek.at(time_int).at(1);
            output[num*num_outputs*max_steps+20+time_int*num_outputs] = time_int;

        }
    }
    py::array ret({number_simulations,max_steps,num_outputs}, output.data() );
    return ret.release();
}



PYBIND11_MODULE(python_sbm, m) {

m.def("run_sbm", &run_sbm, py::return_value_policy::move);
}
