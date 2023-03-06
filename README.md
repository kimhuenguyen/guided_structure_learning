# Code to reproduce the analyses of Nguyen et al. (2023)

In order to run the code in this repository, you need first to install the `learn2count` R package, available at https://github.com/drisso/learn2count.

## Simulation results

To reproduce simulation results, run `scripts/simulation_examples.R`.

Change the structure of the graph to simulate from, by changing the loaded file in line 17. You can control the cardinality, number of simulations, and alpha level by changing the values of lines 24-28.

To see example code of how to create the graph structure that we used for our simulations, see `scripts/samplegraphs.R`.
