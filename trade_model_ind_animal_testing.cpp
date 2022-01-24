#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <iomanip>

using namespace std;

//Structure containing information of individual farms
struct farm {
    int position_in_rate_vectors;

    double eta_equi; //demand generation rate at price equilibrium (read in from file)
    double zeta_equi; //supply generation rate at price equilibrium (read in from file)
    double eta; //current demand generation rate
    double zeta; //current supply generation rate
    double demand = 0.0; //current demand
    double supply = 0.0; //current supply

    vector<int> traders; //list of current trade partnerships
    vector<int> outgoing; //list of farms current selling to

    int disease_state = 0; //current disease state: 0-Susceptible, 1-Infected
    int first_inf = 0;

    double tot_add_rate = 0.0; //current total partnership formation rate for farm
    double tot_delete_rate = 0.0; //current total partnership cessation rate for farm
    double tot_trade_rate = 0.0; //current total trade rate for farm

    double a = 0.0; //partnership formation rate constant (read in from file)
    double d = 0.0; //partnership cessation rate (read in from file)
    double b = 0.0; //trade rate constant (read in from file)

    double base_a = 0.0; //stores value of a before any changes are made via epsilon_a
    double base_d = 0.0; //stores value of d before any changes are made via epsilon_d
    double base_b = 0.0; //stores value of b before any changes are made via epsilon_b

    double sum_traders_equi = 0.0; //stores cumulative number of trading partners after t=50
    double sum_traders_equi_change = 0.0; //stores number of partnership formations after t=50
    double sum_trades_equi = 0.0; //stores cumulative number of trades after t = 50
};

vector<farm> farms; //each object of structure farm is stored as an element of a vector. Each element represents a farm

int main() {
    long double sim_start_time = time(NULL);

    string path(
            "INSERT OUTPUT DIRECTORY HERE");
    string name("INSERT OUTPUT FILE 1 NAME HERE");
    string name2("INSERT OUTPUT FILE 2 NAME HERE");
    string type("INSERT FILE TYPE HERE, e.g. .csv");

    //open file 1 and output column headers
    ofstream out_file;
    out_file.open(path + name + type);
    out_file
            << "time,sensitivity,eta,zeta,needs,supply,traders,trades,trade_size,in_flow,price,infected,income,lost_income"
            << endl;

    //open file 2 and output column headers
    ofstream out_file2;
    out_file2.open(path + name2 + type);
    out_file2 << "sensitivity,eta,zeta,needs,supply,traders,trades,trade_size,in_flow,price,infected,income,lost_income"
              << endl;

    //read in param data
    //READ IN FILE "trade_quantities_final.csv"
    string line;
    ifstream file("INSERT INPUT FILE PATH HERE");

    //check to make sure file is read in correctly
    if (!file.is_open()) {
        cout << "file opening error" << endl;
        exit(EXIT_FAILURE);
    }

    string str;
    getline(file, str);

    //each row of file is separated by ',', each element of row is stored in vector tokens, and assigned to
    //relevant quantity in structure farm
    while (getline(file, line)) {
        istringstream iss{line};

        vector<string> tokens;
        string token;

        while (getline(iss, token, ',')) {
            tokens.push_back(token);
        }

        farm f;
        f.eta_equi = stod(tokens[1]);
        f.eta = stod(tokens[1]);
        f.zeta_equi = stod(tokens[2]);
        f.zeta = stod(tokens[2]);
        f.d = 1 / stod(tokens[4]); //element read in from file is expected partnership duration, hence taking reciprocal
        f.base_d = 1 /
                   stod(tokens[4]); //element read in from file is expected partnership duration, hence taking reciprocal
        f.a = stod(tokens[6]);
        f.base_a = stod(tokens[6]);
        f.b = stod(tokens[7]);
        f.base_b = stod(tokens[7]);
        farms.push_back(f);
    }
    file.close();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    int num_params = 21; //number of parameters to run
    int num_reps = 15; //number of independent reps to run per parameter

    const double N = 1000.0; //number of farms
    const double t_max = 100.0; //length of sim
    const double equi_begin = 25.0; //time to introduce disease
    const double trade_change_time = 50.0; //time at which to consider trading system is in equilibrium and begin
    // calculating equilibrium values for trading quantities

    const double price_equi = 1.0; //price when total supply and demand are equal
    double price; //current price
    double alpha = 0.0001; //price sensitivity parameter
    const double demand_elasticity = 0.412; //price elasticity of demand
    const double supply_elasticity = 0.821; //price elasticity of supply

    const double inf_prob = 0.25; //within-herd prevalence (lambda in main text)
    const double gamma = 1.0 / 3; //infected farm recovery rate
    double B; //probability of infection (function of batch size)

    double m = 0.75; //partnership formation rate exponent (see main text for functional definition of rate)

    double increment_size = 0.1; //size of time step when outputting time series data
    int vector_size = int(t_max / increment_size) + 1; //calculate size of vectors containing time series data

    double sensitivity = 0; //sensitivity of batch test

    for (int np = 0; np < num_params; ++np) {
        long double param_start_time = time(NULL);

        //vectors storing average time series values when averaged over num_reps
        vector<double> avg_avg_eta_unit_time(vector_size, 0);
        vector<double> avg_avg_zeta_unit_time(vector_size, 0);
        vector<double> avg_avg_demand_unit_time(vector_size, 0);
        vector<double> avg_avg_supply_unit_time(vector_size, 0);
        vector<double> avg_avg_traders_unit_time(vector_size, 0);
        vector<double> avg_avg_trades_unit_time(vector_size, 0);
        vector<double> avg_avg_trade_size_unit_time(vector_size, 0);
        vector<double> avg_avg_in_vol_unit_time(vector_size, 0);
        vector<double> avg_avg_price_unit_time(vector_size, 0);
        vector<double> avg_avg_infected_unit_time(vector_size, 0);
        vector<double> avg_avg_income_unit_time(vector_size, 0);
        vector<double> avg_avg_lost_income_unit_time(vector_size, 0);

        sensitivity = 0.05*np;

        for (int nr = 0; nr < num_reps; ++nr) {
            //output consider parameter and rep
            cout << np << '\t' << nr << endl;

            //vectors storing time series data for current rep
            vector<double> avg_eta_unit_time(vector_size, 0);
            vector<double> avg_zeta_unit_time(vector_size, 0);
            vector<double> avg_demand_unit_time(vector_size, 0);
            vector<double> avg_supply_unit_time(vector_size, 0);
            vector<double> avg_traders_unit_time(vector_size, 0);
            vector<double> trade_size_unit_time(vector_size, 0);
            vector<double> trades_unit_time(vector_size, 0);
            vector<double> avg_trade_size_unit_time(vector_size, 0);
            vector<double> avg_trades_unit_time(vector_size, 0);
            vector<double> in_vol_unit_time(vector_size, 0);
            vector<double> avg_in_vol_unit_time(vector_size, 0);
            vector<double> avg_price_unit_time(vector_size, 0);
            vector<double> avg_infected_unit_time(vector_size, 0);
            vector<double> avg_income_unit_time(vector_size, 0);
            vector<double> avg_lost_income_unit_time(vector_size, 0);


            //initialise rng
            unsigned long int seed = time(NULL);
            gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(r, seed);

            //set demand and supply
            //create variables storing system total demand generation rate, supply generation rate, demand, and supply
            double total_eta = 0.0;
            double total_zeta = 0.0;
            double total_demand = 0.0;
            double total_supply = 0.0;

            for (int i = 0; i < N; ++i) {
                total_eta += farms[i].eta_equi;
                total_zeta += farms[i].zeta_equi;
            }
            double avg_eta_init = total_eta / N;
            double avg_zeta_init = total_zeta / N;

            //set initial price to equilibrium price
            price = price_equi;

            //set traders
            //create variable storing total number of edges in system
            double total_traders = 0.0;

            //set disease states
            //system begins disease free. when disease is introduced, only a single farm is infected
            double S = N; //num susceptible
            double I = N - S; //num infected
            int I_init = 1; //num initially infected
            int first_I; //farm index of initial infected

            //set rates
            //create variables storing system total partnership formation rate, cessation rate, trade rate, disease
            //recovery rate, and the sum of these rates
            double total_rate = 0.0;
            double total_add_rate = 0.0;
            double total_delete_rate = 0.0;
            double total_trade_rate = 0.0;
            double total_rec_rate = 0.0;

            //create variables to store changes in rates when an event occurs
            double change_to_demand_rate;
            double change_to_supply_rate;
            double change_to_add_rate;
            double change_to_delete_rate;
            double change_to_trade_rate;

            //here we partition farms into blocks of size sqrt(N). This optimises the number of searches required
            //to calculate what farm performs an event within the Gillespie Algorithm. Rather than checking individual
            //farms, blocks of farms can be checked at the same time.
            int num_partitions = (int) floor(sqrt(N));
            double farms_per_partition = N / num_partitions;

            vector<long double> demand_rates(num_partitions, 0.0);
            vector<long double> supply_rates(num_partitions, 0.0);
            vector<long double> add_rates(num_partitions, 0.0);
            vector<long double> delete_rates(num_partitions, 0.0);
            vector<long double> trade_rates(num_partitions, 0.0);
            vector<int> final_farm(num_partitions, 0);

            for (int i = 0; i < N; ++i) {
                auto position = int(floor(i / farms_per_partition));
                farms[i].position_in_rate_vectors = position;
                demand_rates[position] += farms[i].eta_equi;
                supply_rates[position] += farms[i].zeta_equi;
                final_farm[position] = i;
            }
            total_rate = total_eta + total_zeta; //initially, the total event rate of the system is only the
            //rates at which farms generate supply and demand

            //simulation
            double dt = 0.0; //current time
            double event_rate; //next event rate
            double partial_rate; //value that increments when searching which farm performs event
            int break_check; //variable halting farm searching process once event performing farm is found
            int event_check; //variable that takes value 1 when next event has been found (possibly not needed)

            int demand_check; //variable that takes value 1 when demand generation event occurs
            int supply_check; //variable that takes value 1 when supply generation event occurs
            int add_check; //variable that takes value 1 when partnership formation event occurs
            int delete_check; //variable that takes value 1 when partnership cessation event occurs
            int trade_check; //variable that takes value 1 when trade occurs

            int disease_intro_check = 0; //variable that takes value 1 when disease is introduced
            int inf_check = 0; //variable that takes value 1 if disease dies out before t = 50 (when changes to trade are made).#
            //used to stop and restart sim if disease dies out too early.

            int first_equi_check = 0;

            //gillespie algorithm
            while (dt < t_max) {
                //reset index variables
                event_rate = 0.0;
                partial_rate = 0.0;

                break_check = 0;
                event_check = 0;
                demand_check = 0;
                supply_check = 0;
                add_check = 0;
                delete_check = 0;
                trade_check = 0;

                change_to_demand_rate = 0.0;
                change_to_supply_rate = 0.0;
                change_to_add_rate = 0.0;
                change_to_delete_rate = 0.0;
                change_to_trade_rate = 0.0;

                //first time after equilibrium, store initial equilibrium number of trade partners for each farm
                if (dt >= 50.0 && first_equi_check == 0) {
                    first_equi_check = 1;
                    for (int i = 0; i < N; ++i) {
                        farms[i].sum_traders_equi += farms[i].traders.size();
                        ++farms[i].sum_traders_equi_change;
                    }
                }

                //check to see if disease died out before changes to trade. if disease has, stop sim and restart
                if (dt < 50.0) {
                    if (disease_intro_check == 1) {
                        if (I == 0) {
                            dt = t_max + 1.0;
                            --nr;
                            inf_check = 1;
                            cout << "inf died out" << endl;
                        }
                    }
                }

                //check to see if it's time to introduce infection
                if (disease_intro_check != 1) {
                    if (dt >= equi_begin) {
                        disease_intro_check = 1;
                        while (I < I_init) {
                            auto rand_I = int(gsl_rng_uniform(r) * N); //generate random farm to be first infected
                            if (farms[rand_I].disease_state == 0) {
                                farms[rand_I].disease_state = 1;
                                farms[rand_I].first_inf = 1;
                                ++I; //updated num infected
                                --S; //update num susceptible
                                total_rec_rate = gamma; //update total recovery rate
                                total_rate += gamma; //update total rate
                                first_I = rand_I;

                                //update unit time vectors
                                int time_index =
                                        int(floor(dt / increment_size)) + 1; //calc index of vector given current time
                                if (time_index <= t_max / increment_size) {
                                    avg_infected_unit_time[time_index] += 1.0 / N;
                                }
                            }
                        }
                    }
                }

                //update time
                dt += -log(gsl_rng_uniform_pos(r)) / total_rate;

                //next event rate
                event_rate = gsl_rng_uniform(r) * total_rate;

                //check to see which event occurs
                if (event_rate < total_eta && event_check != 1) {
                    //demand increase
                    event_check = 1;
                    //loop through farm blocks to find which block of farm event occurs in
                    for (int i = 0; i < demand_rates.size() && demand_check != 1; ++i) {
                        partial_rate += demand_rates[i];
                        if (partial_rate >= event_rate) {
                            demand_check = 1;
                            partial_rate -= demand_rates[i];
                            int iterator = 0;
                            if (i > 0) iterator = final_farm[i - 1] + 1; //start next loop from first farm in block
                            //loop through farm block to find which farm performs event
                            for (int it = iterator; it < final_farm[i] + 1 && break_check != 1; ++it) {
                                partial_rate += farms[it].eta;
                                if (partial_rate >= event_rate && break_check != 1) {
                                    break_check = 1;
                                    ++farms[it].demand; //update demand of farm which performs event
                                    ++total_demand; //update total demand of system
                                    double demand_dummy = farms[it].demand;

                                    //update rates
                                    //update trade and partnership formation rates of farm that performs event
                                    for (int j = 0; j < N; ++j) {
                                        if (it != j && farms[j].supply > 0.0) {
                                            double supply_dummy = farms[j].supply;
                                            //check to see if farm j is in list of partners
                                            if (find(farms[it].traders.begin(), farms[it].traders.end(), j) ==
                                                farms[it].traders.end()) {
                                                double old_add_rate = 0.0;
                                                double new_add_rate = 0.0;
                                                old_add_rate = farms[it].a * (demand_dummy - 1) * pow(supply_dummy, m) /
                                                               (N - 1);
                                                new_add_rate =
                                                        farms[it].a * demand_dummy * pow(supply_dummy, m) / (N - 1);
                                                change_to_add_rate += new_add_rate - old_add_rate;
                                                farms[it].tot_add_rate += new_add_rate - old_add_rate;
                                                add_rates[i] += new_add_rate - old_add_rate;
                                            } else {
                                                double old_trade_rate = 0.0;
                                                double new_trade_rate = 0.0;
                                                old_trade_rate = farms[it].b * min(demand_dummy - 1, supply_dummy);
                                                new_trade_rate = farms[it].b * min(demand_dummy, supply_dummy);
                                                change_to_trade_rate += new_trade_rate - old_trade_rate;
                                                farms[it].tot_trade_rate += new_trade_rate - old_trade_rate;
                                                trade_rates[i] += new_trade_rate - old_trade_rate;
                                            }
                                        }
                                    }

                                    //update price
                                    double old_price = price;
                                    price = price_equi * exp(alpha * (total_demand - total_supply));
                                    double price_change = price - old_price;
                                    for (int p = 0; p < N; ++p) {
                                        double old_eta = farms[p].eta;
                                        double new_eta =
                                                farms[p].eta_equi * pow(price / price_equi, -demand_elasticity);
                                        farms[p].eta = new_eta;
                                        change_to_demand_rate += new_eta - old_eta;
                                        demand_rates[farms[p].position_in_rate_vectors] += new_eta - old_eta;
                                        double old_zeta = farms[p].zeta;
                                        double new_zeta =
                                                farms[p].zeta_equi * pow(price / price_equi, supply_elasticity);
                                        farms[p].zeta = new_zeta;
                                        change_to_supply_rate += new_zeta - old_zeta;
                                        supply_rates[farms[p].position_in_rate_vectors] += new_zeta - old_zeta;
                                    }
                                    //update total rates that have been changed by event
                                    total_add_rate += change_to_add_rate;
                                    total_trade_rate += change_to_trade_rate;
                                    total_eta += change_to_demand_rate;
                                    total_zeta += change_to_supply_rate;
                                    total_rate += change_to_add_rate + change_to_trade_rate + change_to_demand_rate +
                                                  change_to_supply_rate;

                                    //update unit time vectors
                                    int time_index = int(floor(dt / increment_size)) + 1;
                                    if (time_index <= t_max / increment_size) {
                                        avg_eta_unit_time[time_index] += change_to_demand_rate / N;
                                        avg_zeta_unit_time[time_index] += change_to_supply_rate / N;
                                        avg_price_unit_time[time_index] += price_change;
                                        avg_demand_unit_time[time_index] += 1.0 / N;
                                    }
                                }
                            }
                        }
                    }
                } else if (event_rate >= total_eta && event_rate < total_eta + total_zeta && event_check != 1) {
                    //supply increase
                    event_check = 1;
                    partial_rate = total_eta;
                    //loop through farm blocks to find which block of farm event occurs in
                    for (int i = 0; i < supply_rates.size() && supply_check != 1; ++i) {
                        partial_rate += supply_rates[i];
                        if (partial_rate >= event_rate) {
                            supply_check = 1;
                            partial_rate -= supply_rates[i];
                            int iterator = 0;
                            if (i > 0) iterator = final_farm[i - 1] + 1; //start next loop from first farm in block
                            //loop through farm block to find which farm performs event
                            for (int it = iterator; it < final_farm[i] + 1 && break_check != 1; ++it) {
                                partial_rate += farms[it].zeta;
                                if (partial_rate >= event_rate && break_check != 1) {
                                    break_check = 1;
                                    ++farms[it].supply; //update supply of farm that performs event
                                    ++total_supply; //update total system supply
                                    double supply_dummy = farms[it].supply;

                                    //update rates
                                    for (int j = 0; j < N; ++j) {
                                        if (it != j && farms[j].demand > 0.0) {
                                            double demand_dummy = farms[j].demand;
                                            //check to see if farm it is in farm j's list of trade partners
                                            if (find(farms[it].outgoing.begin(), farms[it].outgoing.end(), j) ==
                                                farms[it].outgoing.end()) {
                                                double old_add_rate;
                                                double new_add_rate;
                                                old_add_rate =
                                                        farms[j].a * demand_dummy * pow(supply_dummy - 1, m) / (N - 1);
                                                new_add_rate =
                                                        farms[j].a * demand_dummy * pow(supply_dummy, m) / (N - 1);
                                                change_to_add_rate += new_add_rate - old_add_rate;
                                                farms[j].tot_add_rate += new_add_rate - old_add_rate;
                                                add_rates[farms[j].position_in_rate_vectors] +=
                                                        new_add_rate - old_add_rate;
                                            } else {
                                                double old_trade_rate;
                                                double new_trade_rate;
                                                old_trade_rate = farms[j].b * min(demand_dummy, supply_dummy - 1);
                                                new_trade_rate = farms[j].b * min(demand_dummy, supply_dummy);
                                                change_to_trade_rate += new_trade_rate - old_trade_rate;
                                                farms[j].tot_trade_rate += new_trade_rate - old_trade_rate;
                                                trade_rates[farms[j].position_in_rate_vectors] +=
                                                        new_trade_rate - old_trade_rate;
                                            }
                                        }
                                    }

                                    //update price
                                    double old_price = price;
                                    price = price_equi * exp(alpha * (total_demand - total_supply));
                                    double price_change = price - old_price;
                                    for (int p = 0; p < N; ++p) {
                                        double old_eta = farms[p].eta;
                                        double new_eta =
                                                farms[p].eta_equi * pow(price / price_equi, -demand_elasticity);
                                        farms[p].eta = new_eta;
                                        change_to_demand_rate += new_eta - old_eta;
                                        demand_rates[farms[p].position_in_rate_vectors] += new_eta - old_eta;

                                        double old_zeta = farms[p].zeta;
                                        double new_zeta =
                                                farms[p].zeta_equi * pow(price / price_equi, supply_elasticity);
                                        farms[p].zeta = new_zeta;
                                        change_to_supply_rate += new_zeta - old_zeta;
                                        supply_rates[farms[p].position_in_rate_vectors] += new_zeta - old_zeta;
                                    }
                                    //update total rates that have changed due to performed event
                                    total_add_rate += change_to_add_rate;
                                    total_trade_rate += change_to_trade_rate;
                                    total_eta += change_to_demand_rate;
                                    total_zeta += change_to_supply_rate;
                                    total_rate += change_to_add_rate + change_to_trade_rate + change_to_demand_rate +
                                                  change_to_supply_rate;

                                    //update unit time vectors
                                    int time_index = int(floor(dt / increment_size)) + 1;
                                    if (time_index <= t_max / increment_size) {
                                        avg_eta_unit_time[time_index] += change_to_demand_rate / N;
                                        avg_zeta_unit_time[time_index] += change_to_supply_rate / N;
                                        avg_price_unit_time[time_index] += price_change;
                                        avg_supply_unit_time[time_index] += 1.0 / N;
                                    }
                                }
                            }
                        }
                    }
                } else if (event_rate >= total_eta + total_zeta &&
                           event_rate < total_eta + total_zeta + total_add_rate && event_check != 1) {
                    //trader addition
                    event_check = 1;
                    partial_rate = total_eta + total_zeta;
                    //loop through farm blocks to find which block of farm event occurs in
                    for (int i = 0; i < add_rates.size() && add_check != 1; ++i) {
                        partial_rate += add_rates[i];
                        if (partial_rate >= event_rate) {
                            add_check = 1;
                            partial_rate -= add_rates[i];
                            int iterator = 0;
                            if (i > 0) iterator = final_farm[i - 1] + 1; //start next loop from first farm in block
                            //loop through farm block to find which farm performs event
                            for (int it = iterator; it < final_farm[i] + 1 && break_check != 1; ++it) {
                                partial_rate += farms[it].tot_add_rate;
                                if (partial_rate >= event_rate) {
                                    partial_rate -= farms[it].tot_add_rate;
                                    double demand_dummy = farms[it].demand;
                                    for (int j = 0; j < N && break_check != 1; ++j) {
                                        if (it != j && farms[j].supply > 0.0) {
                                            //check to see if farm j is in list of trade partners
                                            if (find(farms[it].traders.begin(), farms[it].traders.end(), j) ==
                                                farms[it].traders.end()) {
                                                partial_rate +=
                                                        farms[it].a * demand_dummy * pow(farms[j].supply, m) / (N - 1);
                                                if (partial_rate >= event_rate && break_check != 1) {
                                                    break_check = 1;
                                                    double supply_dummy = farms[j].supply;
                                                    farms[it].traders.push_back(
                                                            j); //add farm j to list of trade partners
                                                    farms[j].outgoing.push_back(it); //add farm it to list of buyers
                                                    ++total_traders; //update total number of trade partnerships

                                                    //update rates of farm it
                                                    double add_rate;
                                                    double delete_rate;
                                                    double trade_rate;

                                                    add_rate =
                                                            farms[it].a * demand_dummy * pow(supply_dummy, m) / (N - 1);
                                                    delete_rate = farms[it].d;
                                                    trade_rate = farms[it].b * min(demand_dummy, supply_dummy);
                                                    farms[it].tot_add_rate -= add_rate;
                                                    farms[it].tot_delete_rate += delete_rate;
                                                    farms[it].tot_trade_rate += trade_rate;

                                                    add_rates[i] -= add_rate;
                                                    delete_rates[i] += delete_rate;
                                                    trade_rates[i] += trade_rate;

                                                    //update system total rates that have changed due to event
                                                    total_add_rate -= add_rate;
                                                    total_delete_rate += delete_rate;
                                                    total_trade_rate += trade_rate;
                                                    total_rate += trade_rate + delete_rate - add_rate;

                                                    //update unit time vectors
                                                    int time_index = int(floor(dt / increment_size)) + 1;
                                                    if (time_index <= t_max / increment_size) {
                                                        avg_traders_unit_time[time_index] += 1.0 / N;
                                                    }
                                                    //update farm equilibrium values
                                                    if (dt >= 50.0) {
                                                        farms[it].sum_traders_equi += farms[it].traders.size();
                                                        ++farms[it].sum_traders_equi_change;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (event_rate >= total_eta + total_zeta + total_add_rate &&
                           event_rate < total_eta + total_zeta + total_add_rate + total_delete_rate &&
                           event_check != 1) {
                    //trader deletion
                    event_check = 1;
                    partial_rate = total_eta + total_zeta + total_add_rate;
                    //loop through farm blocks to find which block of farm event occurs in
                    for (int i = 0; i < delete_rates.size() && delete_check != 1; ++i) {
                        partial_rate += delete_rates[i];
                        if (partial_rate >= event_rate) {
                            delete_check = 1;
                            partial_rate -= delete_rates[i];
                            int iterator = 0;
                            if (i > 0) iterator = final_farm[i - 1] + 1; //start next loop from first farm in block
                            //loop through farm block to find which farm performs event
                            for (int it = iterator; it < final_farm[i] + 1 && break_check != 1; ++it) {
                                partial_rate += farms[it].tot_delete_rate;
                                if (partial_rate >= event_rate) {
                                    partial_rate -= farms[it].tot_delete_rate;
                                    double demand_dummy = farms[it].demand;
                                    //loop through trade partners of farm it
                                    for (int j = 0; j < farms[it].traders.size() && break_check != 1; ++j) {
                                        double delete_rate;
                                        delete_rate = farms[it].d;
                                        partial_rate += delete_rate;
                                        if (partial_rate >= event_rate && break_check != 1) {
                                            break_check = 1;
                                            int trader_index = farms[it].traders[j];
                                            double supply_dummy = farms[trader_index].supply;

                                            //remove farm j from list of trade partners
                                            farms[it].traders.erase(farms[it].traders.begin() + j);
                                            --total_traders; //update system total number of trade partnerships

                                            //find farm it in farm j's list of buyers and remove
                                            auto position = find(farms[trader_index].outgoing.begin(),
                                                                 farms[trader_index].outgoing.end(), it);
                                            auto position_index = distance(farms[trader_index].outgoing.begin(),
                                                                           position);
                                            farms[trader_index].outgoing.erase(
                                                    farms[trader_index].outgoing.begin() + position_index);

                                            //update rates of farm it
                                            double add_rate;
                                            double trade_rate;
                                            add_rate = farms[it].a * demand_dummy * pow(supply_dummy, m) / (N - 1);
                                            trade_rate = farms[it].b * min(demand_dummy, supply_dummy);
                                            farms[it].tot_delete_rate -= delete_rate;
                                            farms[it].tot_add_rate += add_rate;
                                            farms[it].tot_trade_rate -= trade_rate;

                                            delete_rates[i] -= delete_rate;
                                            add_rates[i] += add_rate;
                                            trade_rates[i] -= trade_rate;

                                            //update system total rates changed by event performed
                                            total_delete_rate -= delete_rate;
                                            total_add_rate += add_rate;
                                            total_trade_rate -= trade_rate;
                                            total_rate += add_rate - trade_rate - delete_rate;

                                            //update unit time vectors
                                            auto time_index = int(floor(dt / increment_size)) + 1;
                                            if (time_index <= t_max / increment_size) {
                                                avg_traders_unit_time[time_index] -= 1.0 / N;
                                            }
                                            //update equilibrium values
                                            if (dt >= 50.0) {
                                                farms[it].sum_traders_equi += farms[it].traders.size();
                                                ++farms[it].sum_traders_equi_change;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (event_rate >= total_eta + total_zeta + total_add_rate + total_delete_rate && event_rate <
                                                                                                        total_eta +
                                                                                                        total_zeta +
                                                                                                        total_add_rate +
                                                                                                        total_delete_rate +
                                                                                                        total_trade_rate &&
                           event_check != 1) {
                    //trade
                    event_check = 1;
                    partial_rate = total_eta + total_zeta + total_add_rate + total_delete_rate;
                    //loop through farm blocks to find which block of farm event occurs in
                    for (int i = 0; i < trade_rates.size() && trade_check != 1; ++i) {
                        partial_rate += trade_rates[i];
                        if (partial_rate >= event_rate) {
                            trade_check = 1;
                            partial_rate -= trade_rates[i];
                            int iterator = 0;
                            if (i > 0) iterator = final_farm[i - 1] + 1; //start next loop from first farm in block
                            //loop through farm block to find which farm performs event
                            for (int it = iterator; it < final_farm[i] + 1 && break_check != 1; ++it) {
                                partial_rate += farms[it].tot_trade_rate;
                                if (partial_rate >= event_rate && break_check != 1) {
                                    partial_rate -= farms[it].tot_trade_rate;
                                    for (int j = 0; j < farms[it].traders.size() && break_check != 1; ++j) {
                                        partial_rate +=
                                                farms[it].b * min(farms[it].demand, farms[farms[it].traders[j]].supply);
                                        if (partial_rate >= event_rate && break_check != 1) {
                                            break_check = 1;

                                            //calc batch size
                                            double trade_size = min(farms[it].demand,
                                                                    farms[farms[it].traders[j]].supply);

                                            double animals_through = trade_size; //this variable will be updated to reflect
                                                                                 //detections of infected animals
                                            double detected_animals = 0.0; //num infected animals detected
                                            double old_price = price;

                                            //update supply of selling farm
                                            total_demand -= trade_size;
                                            total_supply -= trade_size;

                                            //before batch testing is introduced
                                            if (dt < 50) {
                                                farms[it].demand -= trade_size;
                                                total_demand -= trade_size;

                                                //check for infection
                                                if (farms[farms[it].traders[j]].disease_state == 1) {
                                                    if (farms[it].disease_state == 0) {
                                                        B = 1 - pow(1 - inf_prob,
                                                                    trade_size); //calc prob batch infects buyer
                                                        if (gsl_rng_uniform(r) < B) {
                                                            //update num infected and susceptible, and disease state of buyer
                                                            ++I;
                                                            --S;
                                                            farms[it].disease_state = 1;

                                                            //update recovery rates
                                                            total_rec_rate += gamma;
                                                            total_rate += gamma;

                                                            //update unit time vectors
                                                            int time_index = int(floor(dt / increment_size)) + 1;
                                                            if (time_index <= t_max / increment_size) {
                                                                avg_infected_unit_time[time_index] += 1.0 / N;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            //once batch testing is introduced
                                            else{
                                                //check to see if the selling farm is infected. if so, calc number
                                                //of animals in batch that are infected. then, calculate how many
                                                //are detected for infection. animals that test positive are removed,
                                                //animals that do not move on to the buying farm
                                                if (farms[farms[it].traders[j]].disease_state == 1){
                                                    int inf_animals = 0;
                                                    //calc num infected
                                                    for (int ts = 0; ts < trade_size; ++ts){
                                                        if (gsl_rng_uniform(r) < inf_prob){
                                                            ++inf_animals;
                                                        }
                                                    }
                                                    int inf_animals_before = inf_animals;
                                                    //calc num detected
                                                    for (int det = 0; det < inf_animals_before; ++det){
                                                        if (gsl_rng_uniform(r) < sensitivity){
                                                            ++detected_animals;
                                                            --animals_through;
                                                            --inf_animals;
                                                        }
                                                    }
                                                    if (detected_animals > 0.0){
                                                        farms[it].demand -= animals_through;
                                                        total_demand -= animals_through;

                                                        //if infected animals avoid detection, the buying farm becomes
                                                        //infected
                                                        if (inf_animals > 0) {
                                                            if (farms[it].disease_state == 0) {
                                                                farms[it].disease_state = 1;
                                                                //update num infected and susceptible
                                                                ++I;
                                                                --S;
                                                                //update recovery rates
                                                                total_rec_rate += gamma;
                                                                total_rate += gamma;
                                                                //update unit time vectors
                                                                int time_index = int(floor(dt / increment_size)) + 1;
                                                                if (time_index <= t_max / increment_size) {
                                                                    avg_infected_unit_time[time_index] += 1.0 / N;
                                                                }
                                                            }
                                                        }

                                                        //update price - this is the only scenario in which price
                                                        //needs to be updated following a trade as a detection and removal
                                                        //of animals is the only way for trades to be asymmetrical in
                                                        //stock updates of the buying and selling farm
                                                        price = price_equi * exp(alpha * (total_demand - total_supply));
                                                        double price_change = price - old_price;
                                                        for (int p = 0; p < N; ++p) {
                                                            double old_eta = farms[p].eta;
                                                            double new_eta = farms[p].eta_equi *
                                                                             pow(price / price_equi, -demand_elasticity);
                                                            farms[p].eta = new_eta;
                                                            change_to_demand_rate += (new_eta - old_eta);
                                                            demand_rates[farms[p].position_in_rate_vectors] += (new_eta -
                                                                                                               old_eta);

                                                            double old_zeta = farms[p].zeta;
                                                            double new_zeta = farms[p].zeta_equi *
                                                                              pow(price / price_equi,
                                                                                  supply_elasticity);
                                                            farms[p].zeta = new_zeta;
                                                            change_to_supply_rate += (new_zeta - old_zeta);
                                                            supply_rates[farms[p].position_in_rate_vectors] += (
                                                                    new_zeta - old_zeta);
                                                        }
                                                        total_eta += change_to_demand_rate;
                                                        total_zeta += change_to_supply_rate;
                                                        total_rate += change_to_demand_rate + change_to_supply_rate;

                                                        //update unit time vectors
                                                        int time_index = int(floor(dt / increment_size)) + 1;
                                                        if (time_index <= t_max / increment_size) {
                                                            avg_eta_unit_time[time_index] +=
                                                                    change_to_demand_rate / N;
                                                            avg_zeta_unit_time[time_index] +=
                                                                    change_to_supply_rate / N;
                                                            avg_price_unit_time[time_index] += price_change;
                                                        }
                                                    }
                                                    else{
                                                        //no animals are detected
                                                        farms[it].demand -= animals_through;
                                                        total_demand -= animals_through;

                                                        //if there were infected animals in the batch, and if the buying
                                                        //farm is susceptible, the buying farm becomes infected
                                                        if (inf_animals > 0){
                                                            if (farms[it].disease_state == 0){
                                                                farms[it].disease_state = 1;
                                                                //update num infected and susceptible
                                                                ++I;
                                                                --S;
                                                                //update recovery rates
                                                                total_rec_rate += gamma;
                                                                total_rate += gamma;

                                                                //update unit time vectors
                                                                int time_index = int(floor(dt/increment_size)) + 1;
                                                                if (time_index <= t_max / increment_size){
                                                                    avg_infected_unit_time[time_index] += 1.0 / N;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                else{
                                                    //the selling farm is not infected so the whole batch is accepted
                                                    farms[it].demand -= animals_through;
                                                    total_demand -= animals_through;
                                                }
                                            }


                                            //update partnership formation and trade rates relevant to the update of
                                            //demand of farm it
                                            for (int k = 0; k < N; ++k) {
                                                if (k != it) {
                                                    if (farms[k].supply > 0.0) {
                                                        double demand_dummy = farms[it].demand;
                                                        double supply_dummy = farms[k].supply;
                                                        //check to see if farm is a trade partner. if it isn't
                                                        //update formation rate, else update trade rate
                                                        if (find(farms[it].traders.begin(), farms[it].traders.end(),
                                                                 k) == farms[it].traders.end()) {
                                                            double old_add_rate;
                                                            double new_add_rate;
                                                            old_add_rate =
                                                                    farms[it].a * (demand_dummy + animals_through) *
                                                                    pow(supply_dummy, m) / (N - 1);
                                                            new_add_rate =
                                                                    farms[it].a * demand_dummy * pow(supply_dummy, m) /
                                                                    (N - 1);
                                                            change_to_add_rate += new_add_rate - old_add_rate;
                                                            farms[it].tot_add_rate += new_add_rate - old_add_rate;
                                                            add_rates[i] += new_add_rate - old_add_rate;
                                                        } else {
                                                            double old_trade_rate;
                                                            double new_trade_rate;
                                                            if (k == farms[it].traders[j]) {
                                                                old_trade_rate = farms[it].b *
                                                                                 min(demand_dummy + animals_through,
                                                                                     supply_dummy + trade_size);
                                                            } else {
                                                                old_trade_rate = farms[it].b *
                                                                                 min(demand_dummy + animals_through,
                                                                                     supply_dummy);
                                                            }
                                                            new_trade_rate =
                                                                    farms[it].b * min(demand_dummy, supply_dummy);
                                                            change_to_trade_rate += new_trade_rate - old_trade_rate;
                                                            farms[it].tot_trade_rate += new_trade_rate - old_trade_rate;
                                                            trade_rates[i] += new_trade_rate - old_trade_rate;
                                                        }
                                                    }
                                                }
                                            }

                                            //update partnership formation and trade rates relevant to the update of
                                            //supply of selling farm
                                            for (int k = 0; k < N; ++k) {
                                                if (k != farms[it].traders[j]) {
                                                    if (farms[k].demand > 0.0) {
                                                        double demand_dummy = farms[k].demand;
                                                        double supply_dummy = farms[farms[it].traders[j]].supply;
                                                        //check to see if j is a seller. if not, update formation rate,
                                                        //else update trade rate
                                                        if (find(farms[k].traders.begin(), farms[k].traders.end(),
                                                                 farms[it].traders[j]) == farms[k].traders.end()) {
                                                            double old_add_rate;
                                                            double new_add_rate;
                                                            old_add_rate = farms[k].a * demand_dummy *
                                                                           pow(supply_dummy + trade_size, m) / (N - 1);
                                                            new_add_rate =
                                                                    farms[k].a * demand_dummy * pow(supply_dummy, m) /
                                                                    (N - 1);
                                                            change_to_add_rate += new_add_rate - old_add_rate;
                                                            farms[k].tot_add_rate += new_add_rate - old_add_rate;
                                                            add_rates[farms[k].position_in_rate_vectors] +=
                                                                    new_add_rate - old_add_rate;
                                                        } else {
                                                            if (k != it) {
                                                                double old_trade_rate;
                                                                double new_trade_rate;
                                                                old_trade_rate = farms[k].b * min(demand_dummy,
                                                                                                  supply_dummy +
                                                                                                  trade_size);
                                                                new_trade_rate =
                                                                        farms[k].b * min(demand_dummy, supply_dummy);
                                                                change_to_trade_rate += new_trade_rate - old_trade_rate;
                                                                farms[k].tot_trade_rate +=
                                                                        new_trade_rate - old_trade_rate;
                                                                trade_rates[farms[k].position_in_rate_vectors] +=
                                                                        new_trade_rate - old_trade_rate;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            //update system total rates changed by event performed
                                            total_add_rate += change_to_add_rate;
                                            total_trade_rate += change_to_trade_rate;
                                            total_rate += change_to_add_rate + change_to_trade_rate;

                                            //update unit time vectors
                                            auto time_index = int(floor(dt / increment_size)) + 1;
                                            if (time_index <= t_max / increment_size) {
                                                trade_size_unit_time[time_index] += animals_through;
                                                in_vol_unit_time[time_index] += animals_through;
                                                ++trades_unit_time[time_index];
                                                ++avg_trades_unit_time[time_index];

                                                avg_income_unit_time[time_index] += (old_price * animals_through / N);
                                                avg_lost_income_unit_time[time_index] +=
                                                        (old_price * (trade_size - animals_through)) / N;

                                                avg_demand_unit_time[time_index] -= animals_through / N;
                                                avg_supply_unit_time[time_index] -= trade_size / N;
                                            }
                                            //update equilibrium values
                                            if (dt >= 50.0) {
                                                ++farms[it].sum_trades_equi;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (
                        event_rate >= total_eta + total_zeta + total_add_rate + total_delete_rate + total_trade_rate &&
                        event_rate < total_eta + total_zeta + total_add_rate + total_delete_rate + total_trade_rate +
                                     total_rec_rate && event_check != 1) {
                    //recovery
                    event_check = 1;
                    partial_rate = total_eta + total_zeta + total_add_rate + total_delete_rate + total_trade_rate;
                    for (int i = 0; i < N && break_check != 1; ++i) {
                        //check to see if farm is infected
                        if (farms[i].disease_state == 1) partial_rate += gamma;
                        if (partial_rate >= event_rate) {
                            break_check = 1;
                            //update disease state, num infected, and num susceptible
                            farms[i].disease_state = 0;
                            ++S;
                            --I;

                            //update recovery rates
                            total_rec_rate -= gamma;
                            total_rate -= gamma;

                            //update unit time vectors
                            auto time_index = int(floor(dt / increment_size)) + 1;
                            if (time_index <= t_max / increment_size) {
                                avg_infected_unit_time[time_index] -= 1.0 / N;
                            }
                        }
                    }
                }
            }
            //element-wise division to obtain time series average batch sizes
            transform(trade_size_unit_time.begin(), trade_size_unit_time.end(), trades_unit_time.begin(),
                      avg_trade_size_unit_time.begin(), divides<double>());

            //assign values of trade quantities at t=0
            avg_eta_unit_time[0] = avg_eta_init;
            avg_zeta_unit_time[0] = avg_zeta_init;
            avg_price_unit_time[0] = price_equi;
            avg_demand_unit_time[0] = 0.0;
            avg_supply_unit_time[0] = 0.0;
            avg_traders_unit_time[0] = 0.0;
            avg_trade_size_unit_time[0] = 0.0;
            avg_in_vol_unit_time[0] = 0.0;
            avg_infected_unit_time[0] = 0.0;

            if (inf_check != 1) { //only update avg_avg... vectors if disease persists
                for (int i = 0; i < vector_size; ++i) {
                    if (i > 0) {
                        avg_eta_unit_time[i] += avg_eta_unit_time[i - 1];
                        avg_zeta_unit_time[i] += avg_zeta_unit_time[i - 1];
                        avg_demand_unit_time[i] += avg_demand_unit_time[i - 1];
                        avg_supply_unit_time[i] += avg_supply_unit_time[i - 1];
                        avg_traders_unit_time[i] += avg_traders_unit_time[i - 1];
                        avg_price_unit_time[i] += avg_price_unit_time[i - 1];
                        avg_infected_unit_time[i] += avg_infected_unit_time[i - 1];
                        avg_trades_unit_time[i] += avg_trades_unit_time[i - 1];
                        avg_in_vol_unit_time[i] += avg_in_vol_unit_time[i - 1];
                        avg_income_unit_time[i] += avg_income_unit_time[i-1];
                        avg_lost_income_unit_time[i] += avg_lost_income_unit_time[i-1];
                        if (trades_unit_time[i] == 0) avg_trade_size_unit_time[i] = avg_trade_size_unit_time[i - 1];
                    }
                    avg_avg_eta_unit_time[i] += avg_eta_unit_time[i] / num_reps;
                    avg_avg_zeta_unit_time[i] += avg_zeta_unit_time[i] / num_reps;
                    avg_avg_price_unit_time[i] += avg_price_unit_time[i] / num_reps;
                    avg_avg_demand_unit_time[i] += avg_demand_unit_time[i] / num_reps;
                    avg_avg_supply_unit_time[i] += avg_supply_unit_time[i] / num_reps;
                    avg_avg_traders_unit_time[i] += avg_traders_unit_time[i] / num_reps;
                    avg_avg_trade_size_unit_time[i] += avg_trade_size_unit_time[i] / num_reps;
                    avg_avg_in_vol_unit_time[i] += in_vol_unit_time[i] / (N * num_reps);
                    avg_avg_trades_unit_time[i] += trades_unit_time[i] / (N * num_reps);
                    avg_avg_infected_unit_time[i] += avg_infected_unit_time[i] / num_reps;
                    avg_avg_income_unit_time[i] += avg_income_unit_time[i] / num_reps;
                    avg_avg_lost_income_unit_time[i] += avg_lost_income_unit_time[i] / num_reps;
                }
            }

            //reset structure values for next rep
            for (int i = 0; i < N; ++i) {
                farms[i].position_in_rate_vectors = 0;
                farms[i].disease_state = 0;
                farms[i].first_inf = 0;
                farms[i].demand = 0.0;
                farms[i].supply = 0.0;
                farms[i].traders.clear();
                farms[i].outgoing.clear();
                farms[i].tot_add_rate = 0.0;
                farms[i].tot_delete_rate = 0.0;
                farms[i].tot_trade_rate = 0.0;
                farms[i].eta = farms[i].eta_equi;
                farms[i].zeta = farms[i].zeta_equi;
            }
        }

        double avg_eta_equi = 0.0;
        double avg_zeta_equi = 0.0;
        double avg_demand_equi = 0.0;
        double avg_supply_equi = 0.0;
        double avg_traders_equi = 0.0;
        double avg_trades_equi = 0.0;
        double avg_trade_size_equi = 0.0;
        double avg_in_vol_equi = 0.0;
        double avg_price_equi = 0.0;
        double avg_inf_equi = 0.0;
        double avg_income_equi = 0.0;
        double avg_lost_income_equi = 0.0;

        for (int i = 0; i < vector_size; ++i) {
            //output time series data for given epsilon_b value
            out_file << i * increment_size << "," << sensitivity << "," << avg_avg_eta_unit_time[i]
                     << ","
                     << avg_avg_zeta_unit_time[i] << ","
                     << avg_avg_demand_unit_time[i] << "," << avg_avg_supply_unit_time[i] << ","
                     << avg_avg_traders_unit_time[i] << "," << avg_avg_trades_unit_time[i] << ","
                     << avg_avg_trade_size_unit_time[i] << "," << avg_avg_in_vol_unit_time[i] << ","
                     << avg_avg_price_unit_time[i] << ","
                     << avg_avg_infected_unit_time[i] << "," << avg_avg_income_unit_time[i] << ","
                     << avg_avg_lost_income_unit_time[i] << endl;
            //calculate equilibrium average values of trade quantities over final 25 time units of sim
            if (i * increment_size >= t_max - 25.0) {
                avg_eta_equi += avg_avg_eta_unit_time[i] / 25.0;
                avg_zeta_equi += avg_avg_zeta_unit_time[i] / 25.0;
                avg_demand_equi += avg_avg_demand_unit_time[i] / 25.0;
                avg_supply_equi += avg_avg_supply_unit_time[i] / 25.0;
                avg_traders_equi += avg_avg_traders_unit_time[i] / 25.0;
                avg_trades_equi += avg_avg_trades_unit_time[i] / 25.0;
                avg_trade_size_equi += avg_avg_trade_size_unit_time[i] / 25.0;
                avg_in_vol_equi += avg_avg_in_vol_unit_time[i] / 25.0;
                avg_price_equi += avg_avg_price_unit_time[i] / 25.0;
                avg_inf_equi += avg_avg_infected_unit_time[i] / 25.0;
                avg_income_equi += avg_avg_income_unit_time[i] / 25.0;
                avg_lost_income_equi += avg_avg_lost_income_unit_time[i] / 25.0;
            }
        }
        //output equilibrium average trade quantities for given epsilon_b value
        out_file2 << sensitivity << "," << avg_eta_equi << "," << avg_zeta_equi << "," << avg_demand_equi << ","
                  << avg_supply_equi << "," << avg_traders_equi << ","
                  << avg_trades_equi << "," << avg_trade_size_equi << "," << avg_in_vol_equi << "," << avg_price_equi
                  << "," << avg_inf_equi << "," << avg_income_equi << "," << avg_lost_income_equi << endl;

        long double param_end_time = time(NULL);
        cout << "PARAM TIME TAKEN: " << (param_end_time - param_start_time) / 60.0 << " mins" << endl;
    }

    long double sim_end_time = time(NULL);
    cout << "SIM TIME TAKEN: " << (sim_end_time - sim_start_time) / 3600.0 << " hours" << endl;
}
