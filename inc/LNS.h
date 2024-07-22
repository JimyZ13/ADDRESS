#pragma once
#include "BasicLNS.h"
#include "InitLNS.h"
#include <mutex>

//pibt related
#include "simplegrid.h"
#include "pibt_agent.h"
#include "problem.h"
#include "mapf.h"
#include "pibt.h"
#include "pps.h"
#include "winpibt.h"

enum destroy_heuristic { RANDOMAGENTS, RANDOMWALK, INTERSECTION, DESTROY_COUNT, RANDOMWALKPROB};
enum algorithm { CANONICAL, GREEDY, EPSILON, EPSILON_DECAY, UCB, TOPK_GREEDY, TOPK_EPSILON, TOP_EPSILON_DECAY, TOPK_UCB, BERNOULIE, NORMAL };
enum b { BCANONICAL, ADD, REPLACE};
// TODO: adaptively change the neighbor size, that is,
// increase it if no progress is made for a while
// decrease it if replanning fails to find any solutions for several times

class LNS : public BasicLNS
{
public:
    vector<Agent> agents;
    double preprocessing_time = 0;
    double initial_solution_runtime = 0;
    int initial_sum_of_costs = -1;
    int sum_of_costs_lowerbound = -1;
    int sum_of_distances = -1;
    int restart_times = 0;

    LNS(const Instance& instance, double time_limit,
        const string & init_algo_name, const string & replan_algo_name, const string & destroy_name,
        int neighbor_size, int num_of_iterations, bool init_lns, const string & init_destroy_name, bool use_sipp,
        int screen, PIBTPPS_option pipp_option, const string & bandit_algorithm_name, int neighborhoodSizes, int num_agent, string algorithm, double epsilon, double decay, double k, int regions
        , string b);
    ~LNS()
    {
        delete init_lns;
    }
    bool getInitialSolution();
    bool run();
    void validateSolution() const;
    void writeIterStatsToFile(const string & file_name) const;
    void writeResultToFile(const string & file_name) const;
    void writePathsToFile(const string & file_name) const;
    string getSolverName() const override { return "LNS(" + init_algo_name + ";" + replan_algo_name + ")"; }
private:
    InitLNS* init_lns = nullptr;
    string init_algo_name;
    string replan_algo_name;
    bool use_init_lns; // use LNS to find initial solutions
    destroy_heuristic destroy_strategy = RANDOMWALK;
    algorithm algo = GREEDY;
    int num_of_iterations;
    string init_destroy_name;
    PIBTPPS_option pipp_option;
    int num_agent;
    double epsilon;
    double tabu_discount = 0.5; // delayed agent selected probability discount if it's in the tabu_list
    vector<int> delayed_agents;
    vector<int> delay_list;
    b alns_bernoulie = BCANONICAL;
    double decay;
    int k;
    int num_valid_spaces = 0;
    double sub = 0.000025;
    int current_index = -1;
    std::ofstream outFile;
    int regions;
    
    std::vector<int>* q_values;
    std::vector<int>* frequency;

    std::vector<double>* location_q_values;
    std::vector<int>* location_frequency;

    std::vector<int> alpha;
    std::vector<int> beta;
    std::vector<int> location_alpha;
    std::vector<int> location_beta;

    std::vector<double> mu;
    std::vector<double> sigma2;
    std::vector<double> location_mu;
    std::vector<double> location_sigma2;

    PathTable path_table; // 1. stores the paths of all agents in a time-space table;
    // 2. avoid making copies of this variable as much as possible.
    unordered_set<int> tabu_list; // used by randomwalk strategy
    list<int> intersections;

    bool runEECBS();
    bool runCBS();
    bool runPP();
    bool runPIBT();
    bool runPPS();
    bool runWinPIBT();


    MAPF preparePIBTProblem(vector<int>& shuffled_agents);
    void updatePIBTResult(const PIBT_Agents& A, vector<int>& shuffled_agents);

    void chooseDestroyHeuristicbyALNS();

    bool generateNeighborByRandomWalk(int b);
    bool generateNeighborByIntersection();

    int wrapper();
    int location_wrapper();
    int findMostDelayedAgent();
    int greedy();
    int epsilonGreedy();
    int decayEpsilonGreedy();
    int ucb();
    int reverseDecay();
    int epsilonDelay();
    int topKEpsilonGreedy();
    int topKDecayEpsilonGreedy();
    int topKUCB();
    int findRandomAgent() const;
    int location_greedy();
    int location_epsilonGreedy();
    int location_decay_epsilonGreedy();
    int location_topk_epsilon_greedy();
    int location_topk_decay_epsilon_greedy();
    int location_topKUCB();
    int location_UCB();
    int bernoulie();
    int location_bernoulie();
    int normal();
    int location_normal();
    bool generateNeighborByRandomWalkProbSelect();
    int findAgentBasedOnDelay();
    void randomWalk(int agent_id, int start_location, int start_timestep,
                    set<int>& neighbor, int neighbor_size, int upperbound);
};
