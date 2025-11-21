/*
 * Copyright by Paul Rudolph
 * Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
 * https://www.leibniz-hki.de/en/applied-systems-biology.html
 * HKI-Center for Systems Biology of Infection
 * Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
 * Adolf-Reichwein-Straße 23, 07745 Jena, Germany
 *
 * This code is licensed under BSD 2-Clause
 * See the LICENSE file provided with this code for the full license.
 */
#include "external/json.hpp"

#include "utils/io_utils.h"
#include "utils/misc.h"

#include <fstream>
#include <iostream>
#include <utility>

using json = nlohmann::json;
namespace sbm::util {

Parameters get_main_config_parameters(const std::string& main_config) {
    std::ifstream json_file(main_config);
    json json_parameters;
    json_file >> json_parameters;

    // Load main configuration
    Parameters parameters{};
    parameters.general_parameters.main_folder = json_parameters["SBM"]["GeneralParameters"].value("main_folder", "");
    parameters.general_parameters.input_folder = json_parameters["SBM"]["GeneralParameters"].value("input_folder", "");
    parameters.general_parameters.output_folder = json_parameters["SBM"]["GeneralParameters"].value("output_folder", "");

    parameters.general_parameters.standard_seed = json_parameters["SBM"]["GeneralParameters"].value("standard_seed", 42);

    parameters.simulation_parameters.timestep = json_parameters["SBM"]["SimulationParameters"].value("timestep", 1.0);

    parameters.simulation_parameters.timesteps = json_parameters["SBM"]["SimulationParameters"].value("timesteps", 1);
    parameters.model_parameters.resistance_effect = json_parameters["SBM"]["ModelParameters"].value("resistance_effect", false);

    for (const auto& rate : json_parameters["SBM"]["ModelParameters"]["rates"].items()) {
        parameters.model_parameters.rates[rate.key()] = rate.value();
    }

    for (const auto& cell : json_parameters["SBM"]["ModelParameters"]["cells"].items()) {
        parameters.model_parameters.cells[cell.key()] = cell.value();
    }

    json_file.close();
    return parameters;
}
std::vector<double> generate_parameter_vector(sbm::util::ModelParameters model_parameters){
    std::vector<double> internal_parameters(214,0.0);
    internal_parameters.at(0) = round(model_parameters.cells.at("ALIVEEXTERNCANDIDA"));  // 0.0;////00.0;//00.0;//0000.0; //Y,
    internal_parameters.at(1) = round(model_parameters.cells.at("RESISTENTCANDIDA"));    // YR
    internal_parameters.at(2) = round(model_parameters.cells.at("ALIVEIECANDIDA"));      //
    internal_parameters.at(3) = round(model_parameters.cells.at("KILLEDIECANDIDA"));     // 1.0 * pow(10.0, n); //YKE
    internal_parameters.at(4) = round(model_parameters.cells.at("KILLEDEXTERNCANDIDA")); //
    internal_parameters.at(5) = round(model_parameters.cells.at("MONOCYTES"));           // 00.0;//50.0;//4.0;//5833.0;//00000.0; //M00
    internal_parameters.at(105) = round(model_parameters.cells.at("GRANULOCYTES"));      // 247028.0;////00.0;//00.0;//4.0;//3333.0; //G00 = (l+1)*(d+1) + 35.0 * pow(10.0,5);//00.0;//00.0;//4.0;//3333.0; //G00 = (l+1)*(d+1) + 3
    internal_parameters.at(205) = model_parameters.rates.at("PhiG");     // phiG
    internal_parameters.at(206) = model_parameters.rates.at("PhiGStar"); // phiG*
    internal_parameters.at(207) = model_parameters.rates.at("PhiM");     // phiM
    internal_parameters.at(208) = model_parameters.rates.at("KappaM");   // kappaM
    internal_parameters.at(209) = model_parameters.rates.at("KappaG");   // kappaG
    internal_parameters.at(210) = model_parameters.rates.at("Rho");      // rho
    internal_parameters.at(211) = model_parameters.rates.at("Gamma");    // gamma
    internal_parameters.at(212) = model_parameters.rates.at("KappaExt"); // kappaext
    internal_parameters.at(213) = model_parameters.rates.at("GammaR"); // gammaR

    return internal_parameters;
}
} // namespace sbm::util
