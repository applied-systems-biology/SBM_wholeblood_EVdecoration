/*
 * File:   PI_SBM_Solver.h
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

#ifndef PI_SBM_SOLVER_H
#define PI_SBM_SOLVER_H

#include <random>
#include <set>
#include <vector>

#include <utils/io_utils.h>

class PI_SBM_Solver{
  public:
    PI_SBM_Solver(const sbm::util::SimulationParameters& simulation_parameters, std::vector<double>  input_parameter, bool resistance_effect) : model_parameter(std::move(input_parameter)),resistance_effect(resistance_effect) {
        setSimulationProperties(simulation_parameters);
    }
     void set_Parameter(std::vector<double>& parameter);
     void setModel(); // set all model variables that must be initialized

     std::vector<double>& getStartConcentration();
     void getCollectedStates(std::vector<std::vector<double>>& Y, std::vector<std::vector<std::vector<double>>>& M,
                                    std::vector<std::vector<std::vector<double>>>& G, std::vector<std::vector<double>>& R,
                                    std::vector<std::vector<double>>& ICk);

     void simulate(std::minstd_rand& engine);
     void updateYstate(double selected_agent, bool reverse, std::minstd_rand& engine);
     void updateGstate(double selected_agent, bool reverse, std::minstd_rand& engine);
     void updateMstate(double selected_agent, bool reverse, std::minstd_rand& engine);

     int selectYstate(double selected_agent) const;

    bool testForResistance(int Ystate,std::minstd_rand& engine);
    bool testForExtracellularKilling(std::minstd_rand& engine);

     void updateICset();
     void updateICstateVector();
     void updateCurrentYStates();


    static int getSetValueFromIndizes(int i, int j);

     void data_collection();

     void initializeAntimicrobialEffect();
     void initializeResistance();
     void setAntimicrobEffects();
     void setResistance();

    std::vector<std::vector<double>> relativeSubunitEpsilons;

    std::vector<double> start_concentrations;               // Y_temp, G00_temp, M00_temp
    std::vector<std::vector<double>> collected_states;      // all states per timepoint
    std::vector<std::vector<double>> collected_Y_states_no; // includes Y, YR, YKE
    std::vector<std::vector<std::vector<double>>> collected_M_states_no;
    std::vector<std::vector<std::vector<double>>> collected_G_states_no;
    std::vector<std::vector<double>> collected_killedinIC_states_no;

    std::vector<double> number_PhagocytosisEvents;
    std::vector<std::vector<double>> collected_antimicrob_effects;
    std::vector<std::vector<double>> collected_resistance_effects;

    std::vector<std::vector<double>> collected_no_KE; // time, YRA, YRK

    double overall_sum_phag, no_phag_l, no_phag_d, overall_sum_resistance;
    double start_no_M, start_no_G, temp_Y0, start_no_Y_states, M_prob, G_prob, Y_prob; // start_no_Y = temp_Y0
    double temp_Y_states;                                                              // reactive Y states ohne resistance und intracellular gekillten
    double sum_Individuals;
    double declined_phag_events;
    std::vector<std::vector<double>> calculate_MY();
    std::vector<std::vector<double>> calculate_GY();
    std::vector<std::vector<double>> calculate_YL();
    std::vector<std::vector<double>> calculate_YF();
    std::vector<std::vector<double>> calculate_YK();
    std::vector<std::vector<double>> calculate_YKG();
    std::vector<std::vector<double>> calculate_YKM();
    std::vector<std::vector<double>> calculate_YKE();
    std::vector<std::vector<double>> calculate_YAE();
    std::vector<std::vector<double>> calculate_YRK();
    std::vector<std::vector<double>> calculate_YRA();

    std::vector<std::vector<double>> calculate_MMK();
    std::vector<std::vector<double>> calculate_GGK();
    std::vector<std::vector<double>> calculate_MEK();
    std::vector<std::vector<double>> calculate_GEK();

    std::vector<std::vector<double>> calculate_YGA();
    std::vector<std::vector<double>> calculate_YMA();

    std::vector<std::vector<double>> calculate_PE();
    std::vector<std::vector<double>> calculate_kEK();
    std::vector<std::vector<double>> calculate_YR();

  protected:
    bool resistance_effect{};
    std::vector<std::vector<double>> temp_G_states;
    std::vector<std::vector<double>> temp_M_states;
    //    std::vector< std::vector<double> > Y_states;
    double temp_YAE, temp_YKE, temp_YR, temp_YRA, temp_YRK, temp_YAM, temp_YKM, temp_YAG, temp_YKG;
    std::set<int> occupied_G_states;  //(i*10+j) = integer value in set //set with all G states that are occupied with at least one individual
    std::set<int> occupied_M_states;  //(i*10+j) = integer value in set //set with all M states that are occupied with at least one individual
    std::set<int> occupied_GA_states; //(i*10+j) = integer value in set //set with all G states that are occupied with at least one alive individual
    std::set<int> occupied_MA_states; //(i*10+j) = integer value in set //set with all M states that are occupied with at least one alive individual

    double phiG, phiGstar, phiM, kappaM, kappaG, rho, kappaEG, gamma, gammaR;
    double mek, mmk, gek, ggk;

    // 2017-07-16
    bool isEquilibrated;

    // within one timestep
    std::vector<std::vector<double>> prelim_G_states;
    std::vector<std::vector<double>> prelim_M_states;
    //    std::vector< std::vector<double> > Y_states;
    double toAddYKE, toAddYAM, toAddYAG, toAddYKM;
    double toAddYAE, toAddYRK, toAddYRA, toAddYKG;
    double prelim_YAE, prelim_YKE, prelim_YR, prelim_YRA, prelim_YRK, prelim_YAM, prelim_YAG, prelim_YKM, prelim_YKG;
    //    double checked_G, chekced_M;
    std::vector<int> toInsertM;
    std::vector<int> toInsertMA;
    std::vector<int> toInsertG;
    std::vector<int> toInsertGA;

    std::vector<int> toInsertIC;
    std::vector<int> toInsertICA;


    std::pair<int,int> selectIndividualOfICStateM(double selected_agent);
    std::pair<int,int> selectIndividualOfICStateG(double selected_agent);
    bool testForPhagocytosisM(int i, int j, std::minstd_rand& engine);
    bool testForPhagocytosisG(int i, int j, std::minstd_rand& engine);
    bool testForIntracellularKillingM(int i, int j, double probability);
    bool testForIntracellularKillingG(int i, int j, double probability);
    void setSimulationProperties(const sbm::util::SimulationParameters& simulation_parameters);
    double run_time; // complete simulation run time
    double time, time_step;
    double record_data_time_step, record_time;
    int curr_sim_step;
    int no_l, no_d;
    std::vector<double> model_parameter;

    static const bool DEBUG_MODE;
    static const bool REINC_MODE;
    std::string main_directory;


};

#endif /* PI_SBM_SOLVER_H */
