function [] = MasterFunction(n, n_gen, AA_min, N_min, p_p_min, p_c_min, p_t_min, p_d, cat_c, cat_p, w, gamma, gamma2, epsilon, alpha,beta, er, min_length_for_avg,reps)

%addpath('/home/rnunespa/test_job');

tic
%n = Number of cells
%n_gen = Number of time steps
%AA_min = Number of amino acids added at every time step if carbon fixation is uncatalysed 
%N_min = Number of nucleotides added at every time step if carbon fixation is uncatalysed 
%p_p_min = Probability of free polymerisation
%p_c_min = Baseline probability of copying
%p_t_min = Baseline probability of translation
%p_d = Probability of decay
%cat_c = Boolean variable for whether there is catalysis of carbon fixation or not
%cat_p = Boolean variable for whether there is catalysis of copying/translation or not
%w = Boolean variable for whether decay rate is weighed on length or not
%gamma = Catalytic constant that means the minimum length of a catalytic peptide
%epsilon = Catalytic constant meaning the maximum catalytic rate in the length function
%alpha = Catalytic constant for catalysis based on hydrophobicity
%er = error rate in copying and translation 
%min_length_for_avg = minimum length to get average length of polymers from
%reps = Number of repeats of simulations


[n_m_mat_reps, div_per_gen_reps, n_pops_reps, aa_pops_reps, k_avg_vec, avg_length_long_vec, N_m_vec_vec, p_c_vec_vec, meta_data] = selection_function_reps(n, n_gen, AA_min, N_min, p_p_min, p_c_min, p_t_min,p_d,cat_c, cat_p, w,gamma, gamma2, epsilon, alpha, beta, er, min_length_for_avg, reps);



filename = ['results_n_',num2str(n),'_n_gen_',num2str(n_gen),'_AA_min_',num2str(AA_min),'_N_min_',num2str(N_min),'_p_p_min_',num2str(p_p_min),'_p_c_min_',num2str(p_c_min),'_p_t_min_',num2str(p_t_min),'_p_d_',num2str(p_d),'_cat_c_',num2str(cat_c),'_cat_p_',num2str(cat_p),'_w_',num2str(w),'_gamma_',num2str(gamma), '_gamma2_', num2str(gamma2),'_epsilon_',num2str(epsilon),'_alpha_',num2str(alpha),'_beta_',num2str(beta),'_er_',num2str(er),'_min_length_for_avg_',num2str(min_length_for_avg),'_reps_',num2str(reps),'.mat'];

output = {n_m_mat_reps, div_per_gen_reps, n_pops_reps, aa_pops_reps, k_avg_vec, avg_length_long_vec, N_m_vec_vec, p_c_vec_vec, meta_data};

save(filename,'output');
toc
end
