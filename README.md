# Code to reproduce the analyses of Nguyen et al. (2023)

In order to run the code in this repository, you need first to install the `learn2count` R package, available at https://github.com/drisso/learn2count.

## Simulation results

To reproduce simulation results, run `script/simulation_examples.R`.

Change the structure of the graph to simulate from, by changing the loaded file in line 17. You can control the cardinality, number of simulations, and alpha level by changing the values of lines 24-28.

To see example code of how to create the graph structure that we used for our simulations, see `script/samplegraphs.R`.

## Real data

To reproduce results on real data:

1. run `script/prepare_Krasdata.R` to get the data that we considered.
2. run `script/learn_Kras.R` to learn the structure of the graph for the selected subset of the data.

Note that we included the count matrix in the `data/Kras_dataset.RData` object for convenience.
The full, raw data are available at `data/Kras`.
