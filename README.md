# solveOLG closed economy w/o income effects in R (GFT)
Solves a simple AK-OLG-model for a closed economy without income effects in R (using GFT algorithm)

## About
Shows how to solve a simple deterministic overlapping-generations model (OLG) of Auerbauch-Kotlikoff type, solving for the transition path between two steady-states. The code implements the Generalized Fair-Taylor algorithm (Wilcoxen, P. 1990) that iterates over the expectations of agents' until they are aligned with actual outcomes.

A model description can be found here: <https://github.com/solveCGE/solveOLG_doc>. Abstracting from income effects of labor supply allows for an explicit solution of the consumption function.

## How to run
Parameters can be set in `calib.R`. Policy shocks are defined in `run.R` (or just uncomment some of the predefined exemplary shocks). The model is then solved by just running `run.R`. 

## Author
Philip Schuster
