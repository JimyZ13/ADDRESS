# ADDRESS: Bandit-based Adaptive Anytime Multi-Agent Path Finding

This software is based on the latest [MAPF-LNS implementation](https://github.com/Jiaoyang-Li/MAPF-LNS2) from [1] and [Anytime Multi-Agent Path Finding](https://github.com/thomyphan/anytime-mapf) from [2].

ADDRESS is a bandit-enhanced anytime multi-agent path-finding algorithm. The algorithmic framework is based on MAPF-LNS [1] with the neighborhood destroy heuristic replaced by a Thomson Sampling Bandit. 

## Usage
The code requires the external libraries [`BOOST 1.81.0`](https://www.boost.org/) and [`Eigen 3.3`](https://eigen.tuxfamily.org/), and [`CMake`](https://cmake.org) for building the code. 
    
After installing both libraries go to the root folder of this repository and run the following commands: 
```shell script
cmake -DCMAKE_BUILD_TYPE=RELEASE .
make
```

Run the code with:
```
./address -m warehouse-20-40-10-2-2.map -a warehouse-20-40-10-2-2-random-12.scen -o test -k 1000 -t 60 --stats "stats.txt" --outputPaths=paths.txt --seed=0 --maxIterations=1000000 --destroyStrategy=Prob --screen=2 --algorithm=bernoulie --k=32

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
./balance --help
```

We provide example instance files `random-32-32-20.map` and `random-32-32-20-random-1.scen` in the repo. More instances can be downloaded from the [MovingAI MAPF benchmark](https://movingai.com/benchmarks/mapf/index.html).

## Credits

The software is mainly based on code developed by Jiaoyang Li and Zhe Chen in [MAPF-LNS2](https://github.com/Jiaoyang-Li/MAPF-LNS2).

The rule-based MAPF solvers (i.e., PPS, PIBT, and winPIBT) inside the software were borrowed from 
https://github.com/Kei18/pibt/tree/v1.3

BALANCE is released under USC–Research License. See license.txt for further details.

## References

- [1] J. Li et al. *"MAPF-LNS2: Fast Repairing for Multi-Agent Path Finding via Large Neighborhood Search"*. AAAI 2022.
- [2] T. Phan et al. *"Adaptive Anytime Multi-Agent Path Finding using Bandit-Based Large Neighborhood Search"*. AAAI 2024.
