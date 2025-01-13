# An Improved Method for Predicting Local Stability of Microbial Communities with Multiple Interactions - Code Implementation

## Overview

This code package implements the experimental methods presented in the paper *An improved method for predicting local stability of microbial communities with multiple interactions*. The paper introduces a novel approach for predicting the local stability of microbial communities under varying parameter settings. Using this code, users can reproduce the main experimental results from the paper, including the prediction and analysis of microbial community stability. The implementation provides the complete code for the proposed method, supporting comprehensive exploration of community stability analysis.

## Core Method Code

- **funx.m**: Implements the gLV (Lotka-Volterra) dynamical equations.
- **x_ode45.m**: Uses the numerical solver `ode45` to simulate the dynamics of species abundances in the gLV model.
- **BA.m** and **ER.m**: Provide methods for constructing empirical networks (Barabási-Albert and Erdős-Rényi models).
- **DR.m**: Implements the dimensionality reduction algorithm for the dynamics.
- **Tb.m**: Computes the distribution range of the bulk eigenvalues of the Jacobian matrix.
- **C0300.mat, C0301.mat, C0306.mat, C0308.mat**: Real network data used for experimental validation.

## Experiment Code

- **xi_comp.m**: Validates the effectiveness of the dimensionality reduction algorithm on a single node.
- **xm_comp_ER.m**: Validates the effectiveness of the dimensionality reduction algorithm on the entire system.
- **Eig_dist_BA.m**, **Eig_dist_ER.m**: Constructs the Jacobian matrix based on the gLV model and compares the simulated eigenvalue distribution (including λ_bulk and λ_out) with the analytical solution. Users can adjust parameter settings according to their needs, but the parameters must be reasonable.
- **Eig_comp_BA.m**, **Eig_comp_ER.m**: Figure 3b — Validates the effectiveness of the stability analytical methods based on different network models.
- **Eig_xi_BA.m**, **Eig_xi_ER.m**: Figure 4 — Analyzes the relationship between species abundance and community stability. It explores the connection between species abundance and local stability of the community by calculating the community stability at different species abundance levels.
- **diff_dist_BA.m**, **diff_dist_ER.m**: Figure 5 — Compares the stability errors caused by different approximation methods based on species abundance in different empirical network models.
- **real_net.m**, **fig6.m**: Figure 6 — Processes real network data and validates the effectiveness of the stability assessment method on these data.
- **S_N_p.m**, **S_C_mu.m**, **S_sigma.m**: Plots the three subgraphs in Figure 7 and analyzes the influence of different parameters on the matrix stability of the community.


