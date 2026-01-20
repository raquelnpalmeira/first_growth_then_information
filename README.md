# Protocell evolution model (MATLAB)

This repository contains the MATLAB implementation of the individual-based protocell model described in *First Growth, then information: the path to genetic heredity in protocells*.  


---


### Running simulations

Simulations are launched via:

- **`testMatlabFunction.m`**  
  Top-level entry point. Runs a full simulation for a given parameter set and number of stochastic replicates, saving all outputs to disk.

- **`selection_function_reps.m`**  
  Repeats the same parameter set across multiple independent stochastic replicates.

- **`selection_function_test.m`**  
  Runs a single simulation trajectory for a population of protocells, handling time stepping, data recording, and protocell division.

---

### Within-protocell dynamics

All molecular processes acting *within a single protocell* during one time step are implemented in:

- **`model_fun.m`**

This function applies, in sequence: monomer addition, nucleotide polymerisation, RNA copying, RNA translation, and polymer decay.

The following functions implement individual molecular processes used by `model_fun`:

- **Polymerisation**
  - `pairs_to_pol_fun.m`
  - `polfun.m`

- **Copying and translation**
  - `copy_fun.m`
  - `translation_fun.m`
  - `matrix2fvector.m`
  - `intersectWithRepetitions.m`
  - `rank_to_index_fun.m`

- **Decay**
  - `sample_pol_fun.m`
  - `pol_decay_fun.m`

- **Bookkeeping**
  - `avg_length_fun.m`
  - `total_mon_fun.m`

---

### Catalysis

- **`pop_cat_fun.m`**  
  Computes catalytic contributions from peptide populations to CO₂ fixation and templated polymerisation.

---

### Growth, division, and selection

- **`selection_function_test.m`**  
  Handles protocell growth, division thresholds, and replacement dynamics.

- **`div_fun.m`**  
  Implements stochastic partitioning of molecular contents between daughter cells.

---

### Copying and translation errors

- **`error_fun.m`**

The code includes explicit support for copying and translation errors via transition matrices.  
**All simulations reported in the paper used zero error rates (`er = 0`).**

---

## Reproducibility

The model is stochastic; results are reproducible in distribution given identical parameter sets and numbers of replicates.  
Simulations were run in MATLAB (R2021a). Some functions require the Statistics and Machine Learning Toolbox.
