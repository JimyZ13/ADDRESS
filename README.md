# Anytime Multi-Agent Path Finding with an Adaptive Delay-Based Heuristic

This software is based on the latest [MAPF-LNS implementation](https://github.com/Jiaoyang-Li/MAPF-LNS2) from [1] and [Anytime Multi-Agent Path Finding](https://github.com/thomyphan/anytime-mapf) from [2].

Adaptive Delay-based Destroy-and-Repair Enhanced with Success-based Self-Learning (ADDRESS) is a bandit-enhanced anytime multi-agent path-finding algorithm. It is a single-destroy-heuristic variant adaptation of MAPF-LNS proposed in [1] that applies restricted Thompson Sampling to the top-K set of most delayed agents. ADDRESS demonstrates over 50% improvement in AuC on large scale scenarios compared to other state-of-the-art MAPF algorithms. 

## Usage
The code requires the external libraries [`BOOST 1.81.0`](https://www.boost.org/) and [`Eigen 3.3`](https://eigen.tuxfamily.org/), and [`CMake`](https://cmake.org) for building the code. 
    
After installing both libraries go to the root folder of this repository and run the following commands: 
```shell script
cmake -DCMAKE_BUILD_TYPE=RELEASE .
make
```

Run the code with:
```
./address -m Paris_1_256.map -a Paris_1_256-random-1.scen -o test -k 600 -t 60 --stats "stats.txt"
--outputPaths=paths.txt --seed=0 --maxIterations=1000000 --destroyStrategy=RandomWalk --screen=1 --algorithm=bernoulie --k=64 
```

- m: the map file from the MAPF benchmark
- a: the scenario file from the MAPF benchmark
- o: the output file name (no need for file extension)
- k: the number of agents
- t: the runtime limit
- outputPaths: the output file that contains the paths
- algorithm: the algorithm to run for choosing a destroy heuristic in ADDRESS
- b: the flag to insert ADDRESS into MAPF-LNS's ALNS framework 
- epsilon: floating point number to specify the epsilon in epsilon-related algorithms in ADDRESS
- decay: floating point number to specify the decay in decay-related algorithms in ADDRESS
- k: the "k"-most delayed agents which address chooses from in ADDRESS
- regions: the number of regions to split the map into in ADDRESS
- banditAlgo: the multi-armed bandit algorithm used in BALANCE
- neighborCandidateSizes: the number of neighborhood size options

You can find more details and explanations for all parameters with:
```
./address --help
```

The benchmarking test suite for MAPF problems can be downloaded from the [MovingAI MAPF benchmark](https://movingai.com/benchmarks/mapf/index.html).

## Credits

The software is mainly based on code developed by Jiaoyang Li and Zhe Chen in [MAPF-LNS2](https://github.com/Jiaoyang-Li/MAPF-LNS2) and Thomy Phan in [BALANCE](https://github.com/thomyphan/anytime-mapf).

The rule-based MAPF solvers (i.e., PPS, PIBT, and winPIBT) inside the software were borrowed from 
https://github.com/Kei18/pibt/tree/v1.3

## References

- [1] J. Li et al. *"MAPF-LNS2: Fast Repairing for Multi-Agent Path Finding via Large Neighborhood Search"*. AAAI 2022.
- [2] T. Phan et al. *"Adaptive Anytime Multi-Agent Path Finding using Bandit-Based Large Neighborhood Search"*. AAAI 2024.
