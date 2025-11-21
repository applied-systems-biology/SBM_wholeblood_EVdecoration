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

#ifndef SBM_IO_UTILS_H
#define SBM_IO_UTILS_H

#include <external/json.hpp>
#include <array>

using json = nlohmann::json;

namespace sbm::util {

struct GeneralParameters{
    std::string main_folder{};
    std::string input_folder{};
    std::string output_folder{};
    std::string conda_folder{};
    std::uint64_t standard_seed{};
};

struct SimulationParameters{
    double timestep{};
    int timesteps{};
};


struct ModelParameters {
    bool resistance_effect{};
    int number_of_dead_cells{};
    int number_of_parameters{};
    std::map<std::string, double> rates{};
    std::map<std::string, double> cells{};
};


struct Parameters{
    GeneralParameters general_parameters{};
    SimulationParameters simulation_parameters{};
    ModelParameters model_parameters{};
};
Parameters get_main_config_parameters(const std::string& main_config);
std::vector<double> generate_parameter_vector(sbm::util::ModelParameters model_parameters);
} // namespace sbm::util

#endif //SBM_IO_UTILS_H
