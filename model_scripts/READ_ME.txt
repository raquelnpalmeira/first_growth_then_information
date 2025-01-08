MasterFunction takes a set of parameters (what each parameter means is in the comments in the function script), runs a set number of repeats of a simulation, generates and outpur and saves it in a file name containing all the given parameters. 

selection_function_reps runs selection_fun for a set number of repeats. 

selection_function runs model_fun, applying selection in a process similar to a Moran model, for a set number of time steps. 

model_fun contains the dynamics that happen in a single time step. This is, in order, monomer addition, random polymerisation, copying, translation and decay.

