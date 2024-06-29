# ADDRESS: Bandit-based Adaptive Anytime Multi-Agent Path Finding

This software is based on the latest [MAPF-LNS implementation](https://github.com/Jiaoyang-Li/MAPF-LNS2) from [1] and [Anytime Multi-Agent Path Finding](https://github.com/thomyphan/anytime-mapf) from [2].

ADDRESS is a bandit-enhanced anytime multi-agent path-finding algorithm. The algorithmic framework is based on MAPF-LNS [1] and the neighborhood destroy heuristic is replaced with a Thomson Sampling Bandit. 

## Usage
The code requires the external libraries [`BOOST 1.81.0`](https://www.boost.org/) and [`Eigen 3.3`](https://eigen.tuxfamily.org/), and [`CMake`](https://cmake.org) for building the code. 
    
After installing both libraries go to the root folder of this repository and run the following commands: 
```shell script
cmake -DCMAKE_BUILD_TYPE=RELEASE .
make
```

Run the code with:
```
./balance -m Paris_1_256.map -a Paris_1_256-random-1.scen -o test -k 600 -t 60 --stats "stats.txt" --outputPaths=paths.txt --seed=0 --maxIterations=1000000 --destroyStrategy=Adaptive --screen=1 --algorithm=canonical --b=add --epsilon=0.9 --decay=0.005 --k=64 --regions=1000 --banditAlgo=Thompson --neighborCandidateSizes=5

./balance -m random-32-32-20.map -a random-32-32-20-random-1.scen -o test -k 200 -t 60 --stats "stats.txt" --outputPaths=paths.txt --seed=0 --maxIterations=1000000 --destroyStrategy=Adaptive --screen=1 --algorithm=canonical --b=add --epsilon=0.9  --k=32 --regions=1000 --banditAlgo=Roulette --neighborCandidateSizes=1

canonical epsilongreedy=0.3 0.5 0.9 1 bernoulie normal

./balance -m empty-32-32.map -a empty-32-32-random-1.scen -o test -k "300" -t 60  --maxIterations 100000000 --algorithm canonical --screen 1 --epsilon 0 --decay 0  --k 8 --stats "stats.txt" --destroyStrategy "RandomWalk"

./balance -m random-32-32-20.map -a random-32-32-20-random-1.scen -o test -k "200" -t 60  --maxIterations 100000000 --algorithm normal --screen 1 --epsilon 0 --decay 0  --k 8 --stats "stats.txt" --destroyStrategy "Intersection" --regions=64 --neighborSize=8

```

- m: the map file from the MAPF benchmark
- a: the scenario file from the MAPF benchmark
- o: the output file name (no need for file extension)
- k: the number of agents
- t: the runtime limit
- outputPaths: the output file that contains the paths
- banditAlgo: the multi-armed bandit algorithm used in BALANCE
- neighborCandidateSizes: the number of neighborhood size options
- seed: the random seed

-- algorithm: this is the algorithm to use in the destroy strategy 
-- b: this is the control the subset that alns choses from, destroyStrategy must be Adaptive for this to work. This parameter also overrides the algorithm parameter: if b==canonical, algorithm is set to canonical, else algorithm is set to bernoulie. 

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
