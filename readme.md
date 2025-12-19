**Journal Article Materials**

**Title:** *Nest-level phase transition drives synchronized activity bursts in ant colonies*\
**Authors:** Michael Napoli, Simon Garnier, Maurizio Porfiri\
**Submitted To:** *Physical Review X Life*\
**Submitted Date:** December 19, 2025\
**Corresponding Author:** Maurizio Porfiri

---

### File Descriptions

All programs and simulations are run in the Julia programming language. Any file not listed here is a helper file for creating necessary data structures and plots, and for saving data.

* `data/` Contains data necessary to generate figures/perform analysis.
    * `results/` Data created from simulations and saved for figure recreation.
        * `case-0/` Simulations run for constant parameters (Section 3A).
        * `case-1_crit-rho/` Simulations run for variations in the activity transition parameters (Section 3B).
        * `case-1_crit-rho/` Simulations run for variations in the effective temperature and density (Section 3C).

* `gen/` Scripts used to run simulations and save data used for figure/result generation.
    * `figureN.ipynb` File used to generate data used in Figure N.

* `results/` Scripts used to perform analysis and generate the figures included in the manuscript.
    * `figureN.ipynb` File used to construct figure with label N. All data to create figures is saved in the repository, and can be run without running more simulations.

* `model.jl` Helper functions used to simulate the model.

* `network.jl` Helper functions used to create the network between simulated ants.

* `geom.jl` Helper functions used to compute area coverage by simulated ants, pairwise distance, and critical threshold estimates.
