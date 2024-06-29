#include "LNS.h"
#include "ECBS.h"
#include <queue>
#include "LNS.h"
#include "ECBS.h"
#include <queue>
#include <random>
#include <algorithm>
#include <limits>
#include <map>
#include <cmath>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/random_device.hpp>

LNS::LNS(const Instance& instance, double time_limit, const string & init_algo_name, const string & replan_algo_name,
         const string & destroy_name, int neighbor_size, int num_of_iterations, bool use_init_lns,
         const string & init_destroy_name, bool use_sipp, int screen, PIBTPPS_option pipp_option, 
         const string & bandit_algorithm_name, int neighborhoodSizes,
         int num_agent,
         string algorithm, double epsilon, double decay, double k, int regions, string b)  :
         BasicLNS(instance, time_limit, neighbor_size, screen, bandit_algorithm_name, neighborhoodSizes, DESTROY_COUNT),
         init_algo_name(init_algo_name),  replan_algo_name(replan_algo_name), num_of_iterations(num_of_iterations),
         use_init_lns(use_init_lns),init_destroy_name(init_destroy_name),
         path_table(instance.map_size), pipp_option(pipp_option), num_agent(num_agent)
{
    q_values = new std::vector<int>(num_agent, INT_MAX-1);
    frequency = new std::vector<int>(num_agent, 0);
    location_frequency = new std::vector<int>();
    location_q_values = new std::vector<double>();
    start_time = Time::now();
    replan_time_limit = time_limit / 100;
    if (destroy_name == "Adaptive")
    {
        if (b == "canonical"){
            alns_bernoulie = BCANONICAL;
            algorithm = "canonical";
        }
        else if (b == "add"){
            alns_bernoulie = ADD;
        }
        else if (b == "replace"){
            alns_bernoulie = REPLACE;
        }
        int numHeuristics = b == "add" ? DESTROY_COUNT+1 : DESTROY_COUNT;
        ALNS = true;
        heuristicBanditStats.destroy_weights.assign(numHeuristics, 1);
        heuristicBanditStats.destroy_weights_squared.assign(numHeuristics, 1);
        heuristicBanditStats.destroy_counts.assign(numHeuristics, 0);
        decay_factor = 0;
        reaction_factor = 0;
        if(numberOfNeighborhoodSizeCandidates > 0)
        {
            for(int index = 0; index < numHeuristics; index++)
            {
                neighborhoodBanditStats.push_back(new BanditStats());
                neighborhoodBanditStats[index]->destroy_weights.assign(numberOfNeighborhoodSizeCandidates, 1);
                neighborhoodBanditStats[index]->destroy_weights_squared.assign(numberOfNeighborhoodSizeCandidates, 1);
                neighborhoodBanditStats[index]->destroy_counts.assign(numberOfNeighborhoodSizeCandidates, 0);
            }
        }
    }
    else if (destroy_name == "RandomWalk")
        destroy_strategy = RANDOMWALK;
    else if (destroy_name == "Intersection")
        destroy_strategy = INTERSECTION;
    else if (destroy_name == "Random")
        destroy_strategy = RANDOMAGENTS;
    else
    {
        cerr << "Destroy heuristic " << destroy_name << " does not exists. " << endl;
        exit(-1);
    }
     this->epsilon = epsilon;
    this->decay = decay;
    this->k = k;
    if (algorithm == "canonical"){
        algo = CANONICAL;
    }
    else if (algorithm == "greedy"){
        this->epsilon = 0;
        algo = GREEDY;
    } else if (algorithm == "epsilon"){
        algo = EPSILON;
    } else if (algorithm == "decay"){
        algo = EPSILON_DECAY;
    } else if (algorithm == "ucb"){
        algo = UCB;
    } else if (algorithm == "topk-greedy"){
        this->epsilon = 0;
        algo = TOPK_EPSILON;
    } else if (algorithm == "topk-epsilon"){
        algo = TOPK_EPSILON;
    } else if (algorithm == "topk-decay"){
        algo = TOP_EPSILON_DECAY;
    } else if (algorithm == "topk-ucb"){
        algo = TOPK_UCB;
    } else if (algorithm == "bernoulie"){
        algo = BERNOULIE;
    } else if (algorithm == "normal"){
        algo = NORMAL;
    }
    
    this->regions = regions;

    this->neighbor_size = 8;
    alpha = std::vector<int>(num_agent, 1);
    beta = std::vector<int>(num_agent, 1);
    location_alpha = std::vector<int>(regions, 1);
    location_beta = std::vector<int>(regions, 1);
    mu = std::vector<double>(num_agent, 0.0);
    sigma2 = std::vector<double>(num_agent, 1.0);
    location_mu = std::vector<double>();
    location_sigma2 = std::vector<double>();
    int N = instance.getDefaultNumberOfAgents();
    agents.reserve(N);
    for (int i = 0; i < N; i++)
        agents.emplace_back(instance, i, use_sipp);
    preprocessing_time = ((fsec)(Time::now() - start_time)).count();
    if (screen >= 2)
        cout << "Pre-processing time = " << preprocessing_time << " seconds." << endl;
}

bool LNS::run()
{
    // only for statistic analysis, and thus is not included in runtime
    sum_of_distances = 0;
    for (const auto & agent : agents)
    {
        sum_of_distances += agent.path_planner->my_heuristic[agent.path_planner->start_location];
    }

    initial_solution_runtime = 0;
    start_time = Time::now();
    bool succ = getInitialSolution();
    initial_solution_runtime = ((fsec)(Time::now() - start_time)).count();
    if (!succ && initial_solution_runtime < time_limit)
    {
        if (use_init_lns)
        {
            init_lns = new InitLNS(instance, agents, time_limit - initial_solution_runtime,
                    replan_algo_name,init_destroy_name, neighbor_size, screen, bandit_algorithm_name, numberOfNeighborhoodSizeCandidates);
            succ = init_lns->run();
            if (succ) // accept new paths
            {
                path_table.reset();
                for (const auto & agent : agents)
                {
                    path_table.insertPath(agent.id, agent.path);
                }
                init_lns->clear();
                initial_sum_of_costs = init_lns->sum_of_costs;
                sum_of_costs = initial_sum_of_costs;
            }
            initial_solution_runtime = ((fsec)(Time::now() - start_time)).count();
        }
        else // use random restart
        {
            while (!succ && initial_solution_runtime < time_limit)
            {
                succ = getInitialSolution();
                initial_solution_runtime = ((fsec)(Time::now() - start_time)).count();
                restart_times++;
            }
        }
    }

    int searchSuccess = succ? 1 : 0;
    string weights = "";
    iteration_stats.emplace_back(neighbor.agents.size(),
                                 initial_sum_of_costs, initial_solution_runtime, init_algo_name, weights , 0, 0, searchSuccess);
    runtime = initial_solution_runtime;
    if (succ)
    {
        if (screen >= 1)
            cout << "Initial solution cost = " << initial_sum_of_costs << ", "
                 << "runtime = " << initial_solution_runtime << endl;
    }
    else
    {
        cout << "Failed to find an initial solution in "
             << runtime << " seconds after  " << restart_times << " restarts" << endl;
        return false; // terminate because no initial solution is found
    }

    while (runtime < time_limit && iteration_stats.size() <= num_of_iterations)
    {
        runtime =((fsec)(Time::now() - start_time)).count();
        if(screen >= 1)
            validateSolution();
        if (ALNS)
        {
            chooseDestroyHeuristicbyALNS();
        }
        switch (destroy_strategy)
        {
            case RANDOMWALK:
                succ = generateNeighborByRandomWalk(0);
                break;
            case INTERSECTION:
                succ = generateNeighborByIntersection();
                break;
            case RANDOMAGENTS:
                neighbor.agents.resize(agents.size());
                for (int i = 0; i < (int)agents.size(); i++)
                    neighbor.agents[i] = i;
                if (neighbor.agents.size() > neighbor_size)
                {
                    std::random_shuffle(neighbor.agents.begin(), neighbor.agents.end());
                    neighbor.agents.resize(neighbor_size);
                }
                assert(neighbor.agents.size() > 0);
                succ = true;
                break;
            case DESTROY_COUNT:
                succ = generateNeighborByRandomWalk(1);
                break;
            default:
                cerr << "Wrong neighbor generation strategy" << endl;
                exit(-1);
        }
        searchSuccess = succ? 1 : 0;
        if(!succ)
            continue;

        // store the neighbor information
        neighbor.old_paths.resize(neighbor.agents.size());
        neighbor.old_sum_of_costs = 0;
        for (int i = 0; i < (int)neighbor.agents.size(); i++)
        {
            if (replan_algo_name == "PP")
                neighbor.old_paths[i] = agents[neighbor.agents[i]].path;
            path_table.deletePath(neighbor.agents[i], agents[neighbor.agents[i]].path);
            neighbor.old_sum_of_costs += agents[neighbor.agents[i]].path.size() - 1;
        }

        if (replan_algo_name == "EECBS")
            succ = runEECBS();
        else if (replan_algo_name == "CBS")
            succ = runCBS();
        else if (replan_algo_name == "PP")
            succ = runPP();
        else
        {
            cerr << "Wrong replanning strategy" << endl;
            exit(-1);
        }

        if (ALNS) // update destroy heuristics
        {
            const bool condition = neighbor.old_sum_of_costs > neighbor.sum_of_costs;
            double value = (neighbor.old_sum_of_costs - neighbor.sum_of_costs);//
            if(neighbor.agents.size())
            {
                value /= neighbor.agents.size();
            }
            updateDestroyAndNeighborhoodWeights(value, condition);
        }
        runtime = ((fsec)(Time::now() - start_time)).count();
        sum_of_costs += neighbor.sum_of_costs - neighbor.old_sum_of_costs;
        if (screen >= 1)
            cout << "Iteration " << iteration_stats.size() << ", "
                 << "group size = " << neighbor.agents.size() << ", "
                 << "solution cost = " << sum_of_costs << ", "
                 << "remaining time = " << time_limit - runtime << endl;

        double weightSum = 0;
        for (size_t i = 0; i < neighborhoodBanditStats.size(); ++i) {
            for (int j = 0; j < neighborhoodBanditStats[i]->destroy_weights.size(); j++){
                weightSum += neighborhoodBanditStats[i]->destroy_weights[j];
            }
        }
        std::ostringstream oss;
        oss << "\"";
        for (size_t i = 0; i < neighborhoodBanditStats.size(); ++i) {
            string name;
            switch (i)
            {
            case 0:
                name = "random walk";
                break;
            case 1:
                name = "random";
                break;
            case 2:
                name = "intersection";
                break;
            case 3:
                name = "bernoulie";
                break;
            default:
                break;
            }
            oss << name << " : ";
            for (int j = 0; j < neighborhoodBanditStats[i]->destroy_weights.size(); j++){
                oss << j << " = " << neighborhoodBanditStats[i]->destroy_weights[j]/weightSum << "; ";
            }
        }
        oss << "\"";
        std::string weights = oss.str();
        cout << weights<<endl;
        iteration_stats.emplace_back(neighbor.agents.size(), sum_of_costs, runtime, replan_algo_name, weights, 0, 0, searchSuccess);
    }


    average_group_size = - iteration_stats.front().num_of_agents;
    for (const auto& data : iteration_stats)
        average_group_size += data.num_of_agents;
    if (average_group_size > 0)
        average_group_size /= (double)(iteration_stats.size() - 1);

    cout << getSolverName() << ": "
         << "runtime = " << runtime << ", "
         << "iterations = " << iteration_stats.size() << ", "
         << "solution cost = " << sum_of_costs << ", "
         << "initial solution cost = " << initial_sum_of_costs << ", "
         << "failed iterations = " << num_of_failures << endl;
    return true;
}


bool LNS::getInitialSolution()
{
    neighbor.agents.resize(agents.size());
    for (int i = 0; i < (int)agents.size(); i++)
        neighbor.agents[i] = i;
    neighbor.old_sum_of_costs = MAX_COST;
    neighbor.sum_of_costs = 0;
    bool succ = false;
    if (init_algo_name == "EECBS")
        succ = runEECBS();
    else if (init_algo_name == "PP")
        succ = runPP();
    else if (init_algo_name == "PIBT")
        succ = runPIBT();
    else if (init_algo_name == "PPS")
        succ = runPPS();
    else if (init_algo_name == "winPIBT")
        succ = runWinPIBT();
    else if (init_algo_name == "CBS")
        succ = runCBS();
    else
    {
        cerr <<  "Initial MAPF solver " << init_algo_name << " does not exist!" << endl;
        exit(-1);
    }
    if (succ)
    {
        initial_sum_of_costs = neighbor.sum_of_costs;
        sum_of_costs = neighbor.sum_of_costs;
        return true;
    }
    else
    {
        return false;
    }

}

bool LNS::runEECBS()
{
    vector<SingleAgentSolver*> search_engines;
    search_engines.reserve(neighbor.agents.size());
    for (int i : neighbor.agents)
    {
        search_engines.push_back(agents[i].path_planner);
    }

    ECBS ecbs(search_engines, screen - 1, &path_table);
    ecbs.setPrioritizeConflicts(true);
    ecbs.setDisjointSplitting(false);
    ecbs.setBypass(true);
    ecbs.setRectangleReasoning(true);
    ecbs.setCorridorReasoning(true);
    ecbs.setHeuristicType(heuristics_type::WDG, heuristics_type::GLOBAL);
    ecbs.setTargetReasoning(true);
    ecbs.setMutexReasoning(false);
    ecbs.setConflictSelectionRule(conflict_selection::EARLIEST);
    ecbs.setNodeSelectionRule(node_selection::NODE_CONFLICTPAIRS);
    ecbs.setSavingStats(false);
    double w;
    if (iteration_stats.empty())
        w = 5; // initial run
    else
        w = 1.1; // replan
    ecbs.setHighLevelSolver(high_level_solver_type::EES, w);
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime;
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    bool succ = ecbs.solve(T, 0);
    if (succ && ecbs.solution_cost < neighbor.old_sum_of_costs) // accept new paths
    {
        auto id = neighbor.agents.begin();
        for (size_t i = 0; i < neighbor.agents.size(); i++)
        {
            agents[*id].path = *ecbs.paths[i];
            path_table.insertPath(agents[*id].id, agents[*id].path);
            ++id;
        }
        neighbor.sum_of_costs = ecbs.solution_cost;
        if (sum_of_costs_lowerbound < 0)
            sum_of_costs_lowerbound = ecbs.getLowerBound();
    }
    else // stick to old paths
    {
        if (!neighbor.old_paths.empty())
        {
            for (int id : neighbor.agents)
            {
                path_table.insertPath(agents[id].id, agents[id].path);
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        }
        if (!succ)
            num_of_failures++;
    }
    return succ;
}
bool LNS::runCBS()
{
    if (screen >= 2)
        cout << "old sum of costs = " << neighbor.old_sum_of_costs << endl;
    vector<SingleAgentSolver*> search_engines;
    search_engines.reserve(neighbor.agents.size());
    for (int i : neighbor.agents)
    {
        search_engines.push_back(agents[i].path_planner);
    }

    CBS cbs(search_engines, screen - 1, &path_table);
    cbs.setPrioritizeConflicts(true);
    cbs.setDisjointSplitting(false);
    cbs.setBypass(true);
    cbs.setRectangleReasoning(true);
    cbs.setCorridorReasoning(true);
    cbs.setHeuristicType(heuristics_type::WDG, heuristics_type::ZERO);
    cbs.setTargetReasoning(true);
    cbs.setMutexReasoning(false);
    cbs.setConflictSelectionRule(conflict_selection::EARLIEST);
    cbs.setNodeSelectionRule(node_selection::NODE_CONFLICTPAIRS);
    cbs.setSavingStats(false);
    cbs.setHighLevelSolver(high_level_solver_type::ASTAR, 1);
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime; // time limit
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    bool succ = cbs.solve(T, 0);
    if (succ && cbs.solution_cost <= neighbor.old_sum_of_costs) // accept new paths
    {
        auto id = neighbor.agents.begin();
        for (size_t i = 0; i < neighbor.agents.size(); i++)
        {
            agents[*id].path = *cbs.paths[i];
            path_table.insertPath(agents[*id].id, agents[*id].path);
            ++id;
        }
        neighbor.sum_of_costs = cbs.solution_cost;
        if (sum_of_costs_lowerbound < 0)
            sum_of_costs_lowerbound = cbs.getLowerBound();
    }
    else // stick to old paths
    {
        if (!neighbor.old_paths.empty())
        {
            for (int id : neighbor.agents)
            {
                path_table.insertPath(agents[id].id, agents[id].path);
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;

        }
        if (!succ)
            num_of_failures++;
    }
    return succ;
}
bool LNS::runPP()
{
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());
    if (screen >= 2) {
        for (auto id : shuffled_agents)
            cout << id << "(" << agents[id].path_planner->my_heuristic[agents[id].path_planner->start_location] <<
                "->" << agents[id].path.size() - 1 << "), ";
        cout << endl;
    }
    int remaining_agents = (int)shuffled_agents.size();
    auto p = shuffled_agents.begin();
    neighbor.sum_of_costs = 0;
    runtime = ((fsec)(Time::now() - start_time)).count();
    double T = time_limit - runtime; // time limit
    if (!iteration_stats.empty()) // replan
        T = min(T, replan_time_limit);
    auto time = Time::now();

    ConstraintTable constraint_table(instance.num_of_cols, instance.map_size, &path_table);
    while (p != shuffled_agents.end() && ((fsec)(Time::now() - time)).count() < T)
    {
        int id = *p;
        if (screen >= 3)
            cout << "Remaining agents = " << remaining_agents <<
                 ", remaining time = " << T - ((fsec)(Time::now() - time)).count() << " seconds. " << endl
                 << "Agent " << agents[id].id << endl;
        agents[id].path = agents[id].path_planner->findPath(constraint_table);
        if (agents[id].path.empty()) break;
        neighbor.sum_of_costs += (int)agents[id].path.size() - 1;
        if (neighbor.sum_of_costs >= neighbor.old_sum_of_costs)
            break;
        remaining_agents--;
        path_table.insertPath(agents[id].id, agents[id].path);
        ++p;
    }
     if ((destroy_strategy == RANDOMWALK || destroy_strategy == DESTROY_COUNT) && (algo != CANONICAL)){
        if (current_index != -1){
            //update q_values:
            long sum = ((*q_values)[current_index] * ((*frequency)[current_index] - 1)) + sum_of_costs;
            long temp = (sum / (*frequency)[current_index]);
            (*q_values)[current_index] = (int)temp;
            //for bernoulie bandit
            if (algo == BERNOULIE|| destroy_strategy == DESTROY_COUNT){
                int suc = (neighbor.old_sum_of_costs - neighbor.sum_of_costs) > 0 ? 1 : -1;
                if (suc == 1){
                    alpha[current_index]++;
                } else {
                    beta[current_index]++;
                }
            }
            //for normal distribution
            else if (algo == NORMAL){
                double reward = sum_of_costs;
                double n = (*frequency)[current_index];
                double old_mu = mu[current_index];
                double old_sigma2 = sigma2[current_index];

                // Update mean and variance
                double new_mu = old_mu + (reward - old_mu) / n;
                double new_sigma2 = ((n - 1) * old_sigma2 + (reward - old_mu) * (reward - new_mu)) / n;

                mu[current_index] = new_mu;
                sigma2[current_index] = new_sigma2;
            }
        } 
    } else if (destroy_strategy == INTERSECTION){
        if (algo != CANONICAL && current_index != -1){
            //update location q values
            //outFile << "old: " << neighbor.old_sum_of_costs << " new : " << neighbor.sum_of_costs << endl;
            double sum = ((*location_q_values)[current_index] * ((*location_frequency)[current_index] - 1)) + ((neighbor.old_sum_of_costs - neighbor.sum_of_costs));
            double temp = (sum / (*location_frequency)[current_index]);
            (*location_q_values)[current_index] = temp;
            if (algo == BERNOULIE){
                int suc = (neighbor.old_sum_of_costs - neighbor.sum_of_costs) > 0 ? 1 : -1;
                if (suc == 1){
                    location_alpha[current_index]++;
                } else {
                    location_beta[current_index]++;
                }
            }
            else if (algo == NORMAL){
                double reward = sum_of_costs;
                double n = (*location_frequency)[current_index];
                double old_mu = location_mu[current_index];
                double old_sigma2 = location_sigma2[current_index];

                // Update mean and variance
                double new_mu = old_mu + (reward - old_mu) / n;
                double new_sigma2 = ((n - 1) * old_sigma2 + (reward - old_mu) * (reward - new_mu)) / n;

                location_mu[current_index] = new_mu;
                location_sigma2[current_index] = new_sigma2;
            }
        }
    }
    if (remaining_agents == 0 && neighbor.sum_of_costs <= neighbor.old_sum_of_costs) // accept new paths
    {
        return true;
    }
    else // stick to old paths
    {
        if (p != shuffled_agents.end())
            num_of_failures++;
        auto p2 = shuffled_agents.begin();
        while (p2 != p)
        {
            int a = *p2;
            path_table.deletePath(agents[a].id, agents[a].path);
            ++p2;
        }
        if (!neighbor.old_paths.empty())
        {
            p2 = neighbor.agents.begin();
            for (int i = 0; i < (int)neighbor.agents.size(); i++)
            {
                int a = *p2;
                agents[a].path = neighbor.old_paths[i];
                path_table.insertPath(agents[a].id, agents[a].path);
                ++p2;
            }
            neighbor.sum_of_costs = neighbor.old_sum_of_costs;
        }
        return false;
    }
}
bool LNS::runPPS(){
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);
    P.setTimestepLimit(pipp_option.timestepLimit);

    // seed for solver
    auto* MT_S = new std::mt19937(0);
    PPS solver(&P,MT_S);
    solver.setTimeLimit(time_limit);
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}
bool LNS::runPIBT(){
    auto shuffled_agents = neighbor.agents;
     std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);

    // seed for solver
    auto MT_S = new std::mt19937(0);
    PIBT solver(&P,MT_S);
    solver.setTimeLimit(time_limit);
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}
bool LNS::runWinPIBT(){
    auto shuffled_agents = neighbor.agents;
    std::random_shuffle(shuffled_agents.begin(), shuffled_agents.end());

    MAPF P = preparePIBTProblem(shuffled_agents);
    P.setTimestepLimit(pipp_option.timestepLimit);

    // seed for solver
    auto MT_S = new std::mt19937(0);
    winPIBT solver(&P,pipp_option.windowSize,pipp_option.winPIBTSoft,MT_S);
    solver.setTimeLimit(time_limit);
    bool result = solver.solve();
    if (result)
        updatePIBTResult(P.getA(),shuffled_agents);
    return result;
}

MAPF LNS::preparePIBTProblem(vector<int>& shuffled_agents){

    // seed for problem and graph
    auto MT_PG = new std::mt19937(0);
    Graph* G = new SimpleGrid(instance.getMapFile());

    std::vector<Task*> T;
    PIBT_Agents A;

    for (int i : shuffled_agents){
        assert(G->existNode(agents[i].path_planner->start_location));
        assert(G->existNode(agents[i].path_planner->goal_location));
        auto a = new PIBT_Agent(G->getNode( agents[i].path_planner->start_location));
        A.push_back(a);
        Task* tau = new Task(G->getNode( agents[i].path_planner->goal_location));


        T.push_back(tau);
        if(screen>=5){
            cout<<"Agent "<<i<<" start: " <<a->getNode()->getPos()<<" goal: "<<tau->getG().front()->getPos()<<endl;
        }
    }

    return MAPF(G, A, T, MT_PG);

}

void LNS::updatePIBTResult(const PIBT_Agents& A, vector<int>& shuffled_agents){
    int soc = 0;
    for (int i=0; i<A.size();i++){
        int a_id = shuffled_agents[i];

        agents[a_id].path.resize(A[i]->getHist().size());
        int last_goal_visit = 0;
        if(screen>=2)
            std::cout<<A[i]->logStr()<<std::endl;
        for (int n_index = 0; n_index < A[i]->getHist().size(); n_index++){
            auto n = A[i]->getHist()[n_index];
            agents[a_id].path[n_index] = PathEntry(n->v->getId());

            //record the last time agent reach the goal from a non-goal vertex.
            if(agents[a_id].path_planner->goal_location == n->v->getId()
                && n_index - 1>=0
                && agents[a_id].path_planner->goal_location !=  agents[a_id].path[n_index - 1].location)
                last_goal_visit = n_index;

        }
        //resize to last goal visit time
        agents[a_id].path.resize(last_goal_visit + 1);
        if(screen>=2)
            std::cout<<" Length: "<< agents[a_id].path.size() <<std::endl;
        if(screen>=5){
            cout <<"Agent "<<a_id<<":";
            for (auto loc : agents[a_id].path){
                cout <<loc.location<<",";
            }
            cout<<endl;
        }
        path_table.insertPath(agents[a_id].id, agents[a_id].path);
        soc += (int)agents[a_id].path.size()-1;
    }

    neighbor.sum_of_costs =soc;
}

void LNS::chooseDestroyHeuristicbyALNS()
{
    sampleDestroyHeuristicAndNeighborhoodSize();
    if (alns_bernoulie == BCANONICAL ){
        switch (selected_neighbor)
        {
            case 0 : destroy_strategy = RANDOMWALK;cout << "randomwalk" << endl;break;
            case 1 : destroy_strategy = INTERSECTION; cout << "intersection" << endl;break;
            case 2 : destroy_strategy = RANDOMAGENTS; cout << "random" << endl;break;
            default : cerr << "ERROR" << endl; exit(-1);
        }
    }
    else if (alns_bernoulie == REPLACE){
        switch (selected_neighbor)
        {
            case 0 : destroy_strategy = DESTROY_COUNT;cout << "bernoulie" << endl;break;
            case 1 : destroy_strategy = INTERSECTION; cout << "intersection" << endl;break;
            case 2 : destroy_strategy = RANDOMAGENTS; cout << "random" << endl;break;
            default : cerr << "ERROR" << endl; exit(-1);
        }
    }
    else {
        switch (selected_neighbor)
        {
            case 0 : destroy_strategy = RANDOMWALK;cout << "randomwalk" << endl;break;
            case 1 : destroy_strategy = INTERSECTION; cout << "intersection" << endl;break;
            case 2 : destroy_strategy = RANDOMAGENTS; cout << "random" << endl;break;
            case 3 : destroy_strategy = DESTROY_COUNT; cout << "bernoulie" << endl;break;
            default : cerr << "ERROR" << endl; exit(-1);
        }
    }
    
}

bool LNS::generateNeighborByIntersection()
{   
    
    set<int> neighbors_set;
    if (intersections.empty())
    {
        for (int i = 0; i < instance.map_size; i++)
        {
            if (!instance.isObstacle(i) && instance.getDegree(i) > 2)
                intersections.push_back(i);
                num_valid_spaces++;
        }
        for (int i = 0; i < regions; i++){
            (*location_q_values).push_back(INT_MAX-1);
            (*location_frequency).push_back(0); 
            (location_mu).push_back(0.0);
            (location_sigma2).push_back(1.0);
        }
    }
    int location = 0;
    //todo add wrapper and implmeentation
    int region = location_wrapper();
    if (region == -1){
        auto pt = intersections.begin();
        std::advance(pt, rand() % intersections.size());
        location = *pt;
    } else {
        int start = region * (num_valid_spaces / regions);
        int end = start + (num_valid_spaces / regions);
        auto pt = intersections.begin();
        int num_step = 20;
        std::advance(pt, start+(rand() % end));
        location = *pt;
    }
    

    path_table.get_agents(neighbors_set, neighbor_size, location);
    if (neighbors_set.size() < neighbor_size)
    {
        set<int> closed;
        closed.insert(location);
        std::queue<int> open;
        open.push(location);
        while (!open.empty() && (int) neighbors_set.size() < neighbor_size)
        {
            int curr = open.front();
            open.pop();
            for (auto next : instance.getNeighbors(curr))
            {
                if (closed.count(next) > 0)
                    continue;
                open.push(next);
                closed.insert(next);
                if (instance.getDegree(next) >= 3)
                {
                    path_table.get_agents(neighbors_set, neighbor_size, next);
                    if ((int) neighbors_set.size() == neighbor_size)
                        break;
                }
            }
        }
    }
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (neighbor.agents.size() > neighbor_size)
    {
        std::random_shuffle(neighbor.agents.begin(), neighbor.agents.end());
        neighbor.agents.resize(neighbor_size);
    }
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by intersection " << location << endl;
    return true;
}
bool LNS::generateNeighborByRandomWalk(int b)
{
    if (neighbor_size >= (int)agents.size())
    {
        neighbor.agents.resize(agents.size());
        for (int i = 0; i < (int)agents.size(); i++)
            neighbor.agents[i] = i;
        return true;
    }
    int a = -1;
    if (b){
        a = bernoulie();
    }
    else {
        a = wrapper();
    }
    if (a < 0)
        return false;
    
    set<int> neighbors_set;
    neighbors_set.insert(a);
    randomWalk(a, agents[a].path[0].location, 0, neighbors_set, neighbor_size, (int) agents[a].path.size() - 1);
    int count = 0;
    while (neighbors_set.size() < neighbor_size && count < 10)
    {
        int t = rand() % agents[a].path.size();
        randomWalk(a, agents[a].path[t].location, t, neighbors_set, neighbor_size, (int) agents[a].path.size() - 1);
        count++;
        // select the next agent randomly
        int idx = rand() % neighbors_set.size();
        int i = 0;
        for (auto n : neighbors_set)
        {
            if (i == idx)
            {
                a = i;
                break;
            }
            i++;
        }
    }
    if (neighbors_set.size() < 2){
        return false;
    }
    neighbor.agents.assign(neighbors_set.begin(), neighbors_set.end());
    if (screen >= 2)
        cout << "Generate " << neighbor.agents.size() << " neighbors by random walks of agent " << a
             << "(" << agents[a].path_planner->my_heuristic[agents[a].path_planner->start_location]
             << "->" << agents[a].path.size() - 1 << ")" << endl;
    
    return true;
}

int LNS::wrapper(){
    switch(algo){
        case CANONICAL:
            return findMostDelayedAgent();
            break;
        case GREEDY:
            return greedy();
            break;
        case EPSILON:
            return epsilonGreedy();
            break;
        case EPSILON_DECAY:
            return decayEpsilonGreedy();
            break;
        case UCB:
            return ucb();
            break;
        case TOPK_GREEDY:
            return topKEpsilonGreedy();
            break;
        case TOPK_EPSILON:
            return topKEpsilonGreedy();
            break;
        case TOP_EPSILON_DECAY:
            return topKDecayEpsilonGreedy();
            break;
        case TOPK_UCB:
            return topKUCB();
            break;
        case BERNOULIE :
            return bernoulie();
            break;
        case NORMAL :
            return normal();
            break;
    }
}


int LNS::location_wrapper(){
    switch(algo){
        case CANONICAL:
            return -1;
            break;
        case GREEDY:
            return location_greedy();
            break;
        case EPSILON:
            return location_epsilonGreedy();
            break;
        case EPSILON_DECAY:
            return location_decay_epsilonGreedy();
            break;
        case UCB:
            return location_UCB();
            break;
        case BERNOULIE:
            return location_bernoulie();
            break;
        case NORMAL:
            return location_normal();
            break;
        default:
            cout << "Error algorithm not implemented" << endl;
            exit(0);
    }
}



int LNS::location_greedy(){
    double min = INT_MAX;
    int index = 0;
    for (int i =0 ; i < location_q_values->size(); i++){
        if ((*location_q_values)[i] < min ){
            index = i;
            min = (*location_q_values)[i];
        }
    }
    (*location_frequency)[index]++;
    current_index = index;
    return index;
}

int LNS::location_epsilonGreedy(){
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    //if the first selection or epsilon explore
    if (dist(rng) < epsilon) {
        //random arm 
        std::uniform_int_distribution<int> dist2(0, location_q_values->size() - 1);
        int rand = dist2(rng);
        (*location_frequency)[rand]++;
        current_index = rand;
        return rand;
    }
    else{
        //best estimate 
        double max = INT_MIN;
        int index = 0;
        for (int i =0 ; i < location_q_values->size(); i++){
            if ((*location_q_values)[i] > max ){
                index = i;
                max = (*location_q_values)[i];
            }
        }
        (*location_frequency)[index]++;
        current_index = index;
        return index;
    }
}

int LNS::location_decay_epsilonGreedy(){
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    if (epsilon >= 0){
        epsilon -= decay;
    }
    //if the first selection or epsilon explore
    if (dist(rng) < epsilon|| std::all_of(location_frequency->begin(), location_frequency->end(), [](int i){ return i == 0; })) {
        //random arm 
        std::uniform_int_distribution<int> dist2(0, location_q_values->size() - 1);
        int rand = dist2(rng);
        (*location_frequency)[rand]++;
        current_index = rand;
        return rand;
    }
    else{
        //best estimate 
        double max = INT_MIN;
        int index = 0;
        for (int i =0 ; i < location_q_values->size(); i++){
            if ((*location_q_values)[i] >= max ){
                index = i;
                max = (*location_q_values)[i];
            }
        }
        (*location_frequency)[index]++;
        current_index = index;
        return index;
    }
}


int LNS::location_bernoulie(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<double> samples(regions);
    // Sample from the Beta distribution for each agent
    for (int i = 0 ; i < regions ; i++) {
        boost::random::beta_distribution<> beta_dist(location_alpha[i], location_beta[i]);
        samples[i] = beta_dist(gen);
    }

    // Find the agent with the highest sample value
    int index = std::distance(samples.begin(), std::max_element(samples.begin(), samples.end()));
    (*location_frequency)[index]++;
    current_index = index;
    return index;
}

struct LocationComparator{
    bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
        return a.first < b.first;
    }
};

int LNS::findMostDelayedAgent()
{
    int a = -1;
    int max_delays = -1;
    for (int i = 0; i < agents.size(); i++)
    {
        if (tabu_list.find(i) != tabu_list.end())
            continue;
        int delays = agents[i].getNumOfDelays();
        if (max_delays < delays)
        {
            a = i;
            max_delays = delays;
        }
    }
    if (max_delays == 0)
    {
        tabu_list.clear();
        return -1;
    }
    tabu_list.insert(a);
    if (tabu_list.size() == agents.size())
        tabu_list.clear();
    return a;
}

int LNS::greedy()
{
    int max = 0;
    std::mt19937 rng(std::random_device{}());
    //case where all estiamte = 0
    if (std::all_of(frequency->begin(), frequency->end(), [](int i){ return i == 0; })) {
        std::uniform_int_distribution<int> dist2(0, q_values->size() - 1);
        int rand = dist2(rng);
        (*frequency)[rand]++;
        current_index = rand;
        return rand;
    }
    //else choose the best estimate one
    int min = INT_MAX;
    int index = 0;
    for (int i =0 ; i < q_values->size(); i++){
        if ((*q_values)[i] < min ){
            index = i;
            min = (*q_values)[i];
        }
    }
    (*frequency)[index]++;
    current_index = index;
    return index;
}

int LNS::epsilonGreedy()
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    //if the first selection or epsilon explore
        if (dist(rng) < epsilon || std::all_of(frequency->begin(), frequency->end(), [](int i){ return i == 0; })) {
            //random arm 
            std::uniform_int_distribution<int> dist2(0, q_values->size() - 1);
            int rand = dist2(rng);
            (*frequency)[rand]++;
            current_index = rand;
            return rand;
        }
        else{
            //best estimate 
            int min = INT_MAX;
            int index = 0;
            for (int i =0 ; i < q_values->size(); i++){
                if ((*q_values)[i] < min ){
                    index = i;
                    min = (*q_values)[i];
                }
            }
            (*frequency)[index]++;
            current_index = index;
            return index;
        }
}

struct AgentComparator {
    bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
        return a.first < b.first;
    }
};

int LNS::bernoulie() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<double> samples(k);
    //build the top k of most delayed
    
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, AgentComparator> queue;

    // Add all agents to the priority queue
    for (const Agent& agent : agents) {
        int delay = agent.getNumOfDelays();
        int id = agent.id;
        queue.push(std::make_pair(delay, id));
    }

    // Extract the top 10 agents
    std::vector<int> topAgents;
    for (int i = 0; i < k && !queue.empty(); i++) {
        topAgents.push_back(queue.top().second);
        queue.pop();
    }
    // Sample from the Beta distribution for each agent
    for (int i = 0 ; i < k ; i++) {
        boost::random::beta_distribution<> beta_dist(alpha[i], beta[i]);
        samples[i] = beta_dist(gen);
    }

    // Find the agent with the highest sample value
    int best_index = std::distance(samples.begin(), std::max_element(samples.begin(), samples.end()));
    int best_agent = topAgents[best_index];
    (*frequency)[best_agent]++;
    current_index = best_agent;
    return best_agent;
}

int LNS::normal() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<double> samples(k);

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, AgentComparator> queue;

    // Add all agents to the priority queue
    for (const Agent& agent : agents) {
        int delay = agent.getNumOfDelays();
        int id = agent.id;
        queue.push(std::make_pair(delay, id));
    }

    // Extract the top 10 agents
    std::vector<int> topAgents;
    for (int i = 0; i < k && !queue.empty(); i++) {
        topAgents.push_back(queue.top().second);
        queue.pop();
    }

    // Sample from the Gaussian distribution for each agent
    for (int i = 0 ; i < k ; i++) {
        std::normal_distribution<double> normal_dist(mu[i], std::sqrt(sigma2[i]));
        samples[i] = normal_dist(gen);
    }

    // Find the agent with the minimum sample value
    int best_index = std::distance(samples.begin(), std::min_element(samples.begin(), samples.end()));
    int best_agent = topAgents[best_index];
    (*frequency)[best_agent]++;
    current_index = best_agent;
    return best_agent;
}

int LNS::location_normal() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<double> samples(regions);

    // Sample from the Gaussian distribution for each agent
    for (int i = 0 ; i < regions ; i++) {
        std::normal_distribution<double> normal_dist(location_mu[i], std::sqrt(location_sigma2[i]));
        samples[i] = normal_dist(gen);
    }

    // Find the agent with the minimum sample value
    int best_index = std::distance(samples.begin(), std::min_element(samples.begin(), samples.end()));
    (*location_frequency)[best_index]++;
    current_index = best_index;
    return best_index;
}

int LNS::topKEpsilonGreedy()
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    //build the top k of most delayed
    
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, AgentComparator> queue;

    // Add all agents to the priority queue
    for (const Agent& agent : agents) {
        int delay = agent.getNumOfDelays();
        int id = agent.id;
        queue.push(std::make_pair(delay, id));
    }

    // Extract the top 10 agents
    std::vector<int> topAgents;
    for (int i = 0; i < k && !queue.empty(); i++) {
        topAgents.push_back(queue.top().second);
        queue.pop();
    }
    cout << topAgents.size() << endl;
    //if the first selection or epsilon explore
    double i = dist(rng);
    if (dist(rng) < epsilon) {
        //random arm 
        std::uniform_int_distribution<int> dist2(0, topAgents.size() - 1);
        int t = dist2(rng);
        int rand = topAgents[t];
        (*frequency)[rand]++;
        current_index = rand;
        return rand;
    }
    //best estimate among top k agents
    int min = INT_MAX;
    int index = 0;
    for (auto i : topAgents){
        if ((*q_values)[i] < min || (*frequency)[i] == 0 ){
            index = i;
            min = (*q_values)[i];
        }
    }
    (*frequency)[index]++;
    current_index = index;
    return index;
}

int LNS::topKDecayEpsilonGreedy()
{
    std::mt19937 rng(std::random_device{}());
    if (epsilon >= 0){
        epsilon -= decay;
    }
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    //build the top k of most delayed
    
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, AgentComparator> queue;

    // Add all agents to the priority queue
    for (const Agent& agent : agents) {
        int delay = agent.getNumOfDelays();
        int id = agent.id;
        queue.push(std::make_pair(delay, id));
    }

    // Extract the top 10 agents
    std::vector<int> topAgents;
    for (int i = 0; i < k && !queue.empty(); i++) {
        topAgents.push_back(queue.top().second);
        queue.pop();
    }
    //if the first selection or epsilon explore
    if (dist(rng) < epsilon) {
        //random arm 
        std::uniform_int_distribution<int> dist2(0, topAgents.size() - 1);
        int rand = topAgents[dist2(rng)];
        (*frequency)[rand]++;
        current_index = rand;
        return rand;
    }
    //best estimate among top k agents
    int min = INT_MAX;
    int index = 0;
    for (auto i : topAgents){
        if ((*q_values)[i] < min || (*frequency)[i] == 0 ){
            index = i;
            min = (*q_values)[i];
        }
    }
    (*frequency)[index]++;
    current_index = index;
    return index;
}

int LNS::location_UCB()
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    //build the top k of most delayed
    //first explore everything
    for (int i = 0; i < location_frequency->size(); i++){
        if ((*location_frequency)[i] == 0){
            (*location_frequency)[i]++;
            current_index = i;
            return i;
        }
    }
    
    double min = 0;
    int index = 0;
    for (int i = 0 ; i < location_frequency->size(); i++){
        double sqrt = std::sqrt((2*log(10)) / (*location_frequency)[i]);
        double ucb = (*location_q_values)[i] + (*location_q_values)[i]*sqrt;
        if (ucb < min || min == 0){
            index = i;
            min = ucb;
        }
    }
    (*location_frequency)[index]++;
    current_index = index;
    return index;
}



int LNS::topKUCB()
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    //build the top k of most delayed
    //first explore everything
    for (int i = 0; i < frequency->size(); i++){
        if ((*frequency)[i] == 0){
            (*frequency)[i]++;
            current_index = i;
            return i;
        }
    }
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, AgentComparator> queue;

    // Add all agents to the priority queue
    for (const Agent& agent : agents) {
        int delay = agent.getNumOfDelays();
        int id = agent.id;
        queue.push(std::make_pair(delay, id));
    }

    // Extract the top 10 agents
    std::vector<int> topAgents;
    for (int i = 0; i < 10 && !queue.empty(); i++) {
        topAgents.push_back(queue.top().second);
        queue.pop();
    }

    
    double min = 0;
    int index = 0;
    for (auto i : topAgents){
        double sqrt = std::sqrt((2*log(10)) / (*frequency)[i]);
        double ucb = (*q_values)[i] + initial_sum_of_costs*sqrt;
        if (ucb < min || min == 0){
            index = i;
            min = ucb;
        }
    }
    (*frequency)[index]++;
    current_index = index;
    return index;
}

int LNS::decayEpsilonGreedy()
{
    std::mt19937 rng(std::random_device{}());
    cout << "epsilon value : " << epsilon << endl;
    if (epsilon >= 0){
        epsilon -= decay;
    }
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    //if the first selection or epsilon explore
        if (dist(rng) < epsilon || std::all_of(frequency->begin(), frequency->end(), [](int i){ return i == 0; })) {
            //random arm 
            std::uniform_int_distribution<int> dist2(0, q_values->size() - 1);
            int rand = dist2(rng);
            (*frequency)[rand]++;
            current_index = rand;
            return rand;
        }
        else{
            //best estimate 
             int min = INT_MAX;
            int index = 0;
            for (int i =0 ; i < q_values->size(); i++){
                if ((*q_values)[i] < min ){
                    index = i;
                    min = (*q_values)[i];
                }
            }
            (*frequency)[index]++;
            current_index = index;
            return index;
        }
}

int LNS::epsilonDelay(){
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    //if the first selection or epsilon explore
        if (dist(rng) < epsilon || std::all_of(frequency->begin(), frequency->end(), [](int i){ return i == 0; })) {
            //random arm 
            std::uniform_int_distribution<int> dist(0, q_values->size() - 1);
            int rand = dist(rng);
            (*frequency)[rand]++;
            current_index = rand;
            return rand;
        }
        else{
            //best estimate 
            int index = findMostDelayedAgent();
            (*frequency)[index]++;
            current_index = index;
            return index;
        }
}

int LNS::ucb()
{   
    int max = 0;
    //must explore every agent once 
    for (int i = 0; i < frequency->size(); i++){
        if ((*frequency)[i] == 0){
            (*frequency)[i]++;
            current_index = i;
            return i;
        }
    }
    double min = 0;
    int index = 0;
    for (int i =0 ; i < q_values->size(); i++){
        double sqrt = std::sqrt((2*log(10)) / (*frequency)[i]);
        double ucb = (*q_values)[i] + initial_sum_of_costs * sqrt;
        if (ucb < min || min == 0){
            index = i;
            min = ucb;
        }
    }
    (*frequency)[index]++;
    current_index = index;
    return index;

}


int LNS::findRandomAgent() const
{
    int a = 0;
    int pt = rand() % (sum_of_costs - sum_of_distances) + 1;
    int sum = 0;
    for (; a < (int) agents.size(); a++)
    {
        sum += agents[a].getNumOfDelays();
        if (sum >= pt)
            break;
    }
    assert(sum >= pt);
    return a;
}

// a random walk with path that is shorter than upperbound and has conflicting with neighbor_size agents
void LNS::randomWalk(int agent_id, int start_location, int start_timestep,
                     set<int>& conflicting_agents, int neighbor_size, int upperbound)
{
    int loc = start_location;
    for (int t = start_timestep; t < upperbound; t++)
    {
        auto next_locs = instance.getNeighbors(loc);
        next_locs.push_back(loc);
        while (!next_locs.empty())
        {
            int step = rand() % next_locs.size();
            auto it = next_locs.begin();
            advance(it, step);
            int next_h_val = agents[agent_id].path_planner->my_heuristic[*it];
            if (t + 1 + next_h_val < upperbound) // move to this location
            {
                path_table.getConflictingAgents(agent_id, conflicting_agents, loc, *it, t + 1);
                loc = *it;
                break;
            }
            next_locs.erase(it);
        }
        if (next_locs.empty() || conflicting_agents.size() >= neighbor_size)
            break;
    }
}

void LNS::validateSolution() const
{
    int sum = 0;
    for (const auto& a1_ : agents)
    {
        if (a1_.path.empty())
        {
            cerr << "No solution for agent " << a1_.id << endl;
            exit(-1);
        }
        else if (a1_.path_planner->start_location != a1_.path.front().location)
        {
            cerr << "The path of agent " << a1_.id << " starts from location " << a1_.path.front().location
                << ", which is different from its start location " << a1_.path_planner->start_location << endl;
            exit(-1);
        }
        else if (a1_.path_planner->goal_location != a1_.path.back().location)
        {
            cerr << "The path of agent " << a1_.id << " ends at location " << a1_.path.back().location
                 << ", which is different from its goal location " << a1_.path_planner->goal_location << endl;
            exit(-1);
        }
        for (int t = 1; t < (int) a1_.path.size(); t++ )
        {
            if (!instance.validMove(a1_.path[t - 1].location, a1_.path[t].location))
            {
                cerr << "The path of agent " << a1_.id << " jump from "
                     << a1_.path[t - 1].location << " to " << a1_.path[t].location
                     << " between timesteps " << t - 1 << " and " << t << endl;
                exit(-1);
            }
        }
        sum += (int) a1_.path.size() - 1;
        for (const auto  & a2_: agents)
        {
            if (a1_.id >= a2_.id || a2_.path.empty())
                continue;
            const auto & a1 = a1_.path.size() <= a2_.path.size()? a1_ : a2_;
            const auto & a2 = a1_.path.size() <= a2_.path.size()? a2_ : a1_;
            int t = 1;
            for (; t < (int) a1.path.size(); t++)
            {
                if (a1.path[t].location == a2.path[t].location) // vertex conflict
                {
                    cerr << "Find a vertex conflict between agents " << a1.id << " and " << a2.id <<
                            " at location " << a1.path[t].location << " at timestep " << t << endl;
                    exit(-1);
                }
                else if (a1.path[t].location == a2.path[t - 1].location &&
                        a1.path[t - 1].location == a2.path[t].location) // edge conflict
                {
                    cerr << "Find an edge conflict between agents " << a1.id << " and " << a2.id <<
                         " at edge (" << a1.path[t - 1].location << "," << a1.path[t].location <<
                         ") at timestep " << t << endl;
                    exit(-1);
                }
            }
            int target = a1.path.back().location;
            for (; t < (int) a2.path.size(); t++)
            {
                if (a2.path[t].location == target)  // target conflict
                {
                    cerr << "Find a target conflict where agent " << a2.id << " (of length " << a2.path.size() - 1<<
                         ") traverses agent " << a1.id << " (of length " << a1.path.size() - 1<<
                         ")'s target location " << target << " at timestep " << t << endl;
                    exit(-1);
                }
            }
        }
    }
    if (sum_of_costs != sum)
    {
        cerr << "The computed sum of costs " << sum_of_costs <<
             " is different from the sum of the paths in the solution " << sum << endl;
        exit(-1);
    }
}

void LNS::writeIterStatsToFile(const string & file_name) const
{
    if (init_lns != nullptr)
    {
        init_lns->writeIterStatsToFile(file_name + "-initLNS.csv");
    }
    if (iteration_stats.size() <= 1)
        return;
    string name = file_name;
    if (use_init_lns or num_of_iterations > 0)
        name += "-LNS.csv";
    else
        name += "-" + init_algo_name + ".csv";
    std::ofstream output;
    output.open(name);
    // header
    cout << "======"<< ALNS << "=======" << endl;
    if (1 > 0){
        output << "num of agents," <<
            "sum of costs," <<
            "weights," << 
            "runtime," <<
            "cost lowerbound," <<
            "sum of distances," <<
            "MAPF algorithm" << endl;
            for (const auto &data : iteration_stats)
            {
                output << data.num_of_agents << "," <<
                    data.sum_of_costs << "," << data.weights << "," << data.runtime << "," <<
                    max(sum_of_costs_lowerbound, sum_of_distances) << "," <<
                    sum_of_distances << "," <<
                    data.algorithm << endl;
            }
    }
    else {
        output << "num of agents," <<
            "sum of costs," <<
            "runtime," <<
            "cost lowerbound," <<
            "sum of distances," <<
            "MAPF algorithm" << endl;
        for (const auto &data : iteration_stats)
        {
            output << data.num_of_agents << "," <<
                data.sum_of_costs << "," <<
                data.runtime << "," <<
                max(sum_of_costs_lowerbound, sum_of_distances) << "," <<
                sum_of_distances << "," <<
                data.algorithm << endl;
        }
    }
    
    output.close();
}

void LNS::writeResultToFile(const string & file_name) const
{
    if (init_lns != nullptr)
    {
        init_lns->writeResultToFile(file_name + "-initLNS.csv", sum_of_distances, preprocessing_time);
    }
    string name = file_name;
    if (use_init_lns or num_of_iterations > 0)
        name += "-LNS.csv";
    else
        name += "-" + init_algo_name + ".csv";
    std::ifstream infile(name);
    bool exist = infile.good();
    infile.close();
    if (!exist)
    {
        ofstream addHeads(name);
        addHeads << "runtime,solution cost,initial solution cost,lower bound,sum of distance," <<
                 "iterations," <<
                 "group size," <<
                 "runtime of initial solution,restart times,area under curve," <<
                 "LL expanded nodes,LL generated,LL reopened,LL runs," <<
                 "preprocessing runtime,solver name,instance name,success,selected_neighbor,neighbor_size" << endl;
        addHeads.close();
    }
    uint64_t num_LL_expanded = 0, num_LL_generated = 0, num_LL_reopened = 0, num_LL_runs = 0;
    for (auto & agent : agents)
    {
        agent.path_planner->reset();
        num_LL_expanded += agent.path_planner->accumulated_num_expanded;
        num_LL_generated += agent.path_planner->accumulated_num_generated;
        num_LL_reopened += agent.path_planner->accumulated_num_reopened;
        num_LL_runs += agent.path_planner->num_runs;
    }
    double auc = 0;
    if (!iteration_stats.empty())
    {
        auto prev = iteration_stats.begin();
        auto curr = prev;
        ++curr;
        while (curr != iteration_stats.end() && curr->runtime < time_limit)
        {
            auc += (prev->sum_of_costs - sum_of_distances) * (curr->runtime - prev->runtime);
            prev = curr;
            ++curr;
        }
        auc += (prev->sum_of_costs - sum_of_distances) * (time_limit - prev->runtime);
    }
    ofstream stats(name, std::ios::app);
    stats << runtime << "," << sum_of_costs << "," << initial_sum_of_costs << "," <<
          sum_of_costs_lowerbound << "," << sum_of_distances << "," <<
          iteration_stats.size() << "," << average_group_size << "," <<
          initial_solution_runtime << "," << restart_times << "," << auc << "," <<
          num_LL_expanded << "," << num_LL_generated << "," << num_LL_reopened << "," << num_LL_runs << "," <<
          preprocessing_time << "," << getSolverName() << "," << instance.getInstanceName() << "," << iteration_stats.back().success << "," << selected_neighbor << "," << neighbor_size << endl;
    stats.close();
}

void LNS::writePathsToFile(const string & file_name) const
{
    std::ofstream output;
    output.open(file_name);

    for (const auto &agent : agents)
    {
        output << "Agent " << agent.id << ":";
        for (const auto &state : agent.path)
            output << "(" << instance.getRowCoordinate(state.location) << "," <<
                            instance.getColCoordinate(state.location) << ")->";
        output << endl;
    }
    output.close();
}
