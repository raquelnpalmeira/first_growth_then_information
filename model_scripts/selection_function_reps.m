function [n_m_mat_reps, div_per_gen_reps, n_pops_reps, aa_pops_reps, k_avg_vec, avg_length_long_vec, N_m_vec_vec, p_c_vec_vec, meta_data] = selection_function_reps(n, n_gen, AA_min, N_min, p_p_min, p_c_min, p_t_min,p_d,cat_c, cat_p, w,gamma, gamma2, epsilon, alpha,beta,er, min_length_for_avg, reps)

n_m_mat_reps = cell(1,reps);
div_per_gen_reps = cell(1,reps);
n_pops_reps = cell(1,reps);
aa_pops_reps = cell(1,reps);
k_avg_vec = cell(1,reps);
avg_length_long_vec = cell(1,reps);
N_m_vec_vec = cell(1,reps);
p_c_vec_vec = cell(1,reps);

% Number of workers in the parallel pool
n_w = 4; 

% Create parallel pool
parpool(n_w); 

% Set maximum number of workers to be used by the parfor loop (use one less than the maximum so that matlab overhead is still run by one core) 
M = n_w - 1;


parfor (i = 1:reps, M)
%for i = 1:reps
[n_m_mat_reps{1,i}, div_per_gen_reps{1,i}, n_pops_reps{1,i}, aa_pops_reps{1,i}, k_avg_vec{1,i}, avg_length_long_vec{1,i}, N_m_vec_vec{1,i}, p_c_vec_vec{1,i}] = selection_function(n, n_gen, AA_min, N_min, p_p_min, p_c_min, p_t_min,p_d,cat_c, cat_p, w,gamma, gamma2, epsilon, alpha, beta,er, min_length_for_avg);
end 


Parameter = ["n"; "n_gen"; "AA_min"; "N_min"; "p_p_min"; "p_c_min"; "p_t_min";"p_d";"cat_c"; "cat_p"; "w";"gamma";"gamma2"; "epsilon"; "alpha";"beta";"er"; "min_length_for_avg"; "reps"];
Value = [n, n_gen, AA_min, N_min, p_p_min, p_c_min, p_t_min,p_d,cat_c, cat_p, w,gamma,gamma2, epsilon, alpha,beta,er, min_length_for_avg, reps]';
meta_data = table(Parameter, Value);
%meta_data = table([n, n_gen, AA_min, N_min, p_p_min, p_c_min, p_t_min,p_d,cat_c, cat_p, w,gamma, epsilon, alpha,er, min_length_for_avg, reps]);

end
