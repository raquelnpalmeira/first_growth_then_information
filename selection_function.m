% Function for evolutionary model
function [n_m_mat, div_per_gen, n_pops, aa_pops, k_avg_vec, avg_length_long_vec, N_m_vec_vec, p_c_vec_vec] = selection_function(n, n_gen, AA_min, N_min, p_p_min, p_c_min, p_t_min,p_d,cat_c, cat_p, w,gamma, gamma2, epsilon, alpha, beta, er, min_length_for_avg)
%% Run x cells in parallel
t = 1; % number of time steps per generation
%n_gen = 1000; 
%n  = 50; % population size (use an even number)
n_rec = 5;
size_rec_int = n_gen/n_rec;
threshold_size = 1e4;%1000; %size at which division happens
max_l = 40;%40; % maximum length of polymers (dimentions for the population of molecules)
n1 = zeros(max_l,max_l+1);
l_vec = 1:length(n1(:,1));


%initialise population cell arrays
aa_pop = cell(1,n);
n_pop = cell(1,n);


%fill initial population cell arrays with empty cells

for i = 1:n
    aa_pop{i} = n1;
    n_pop{i} = n1;
end

%% set up parameters

% catalytic constants
%gamma = 6;%12; % 
% epsilon = 1.5;%2e-3;%1e-3;%2e-3;%%3;%[0,3,6,12];%1.5; % 
% alpha = 4;%0.8;% try higher values so dimers better than two monomers % catalytic constant for catalysis based on length



% p_p_min = 0.02;% 0.02;%0.03;%[0,0.03,0.06,0.09];%0.05;
% p_c_min = 0.1;%0.05;%0.02;%0;%0.005;%0.005;%0.05;
% p_t_min = 0.1;%0.01;%0.005;%0.005;%[0,0.001,0.01,0.1];%0;%
%p_d = 0.001;%[0,0.001,0.01,0.1];%0;%0.005;%
% AA_min = 1e4;%5e3;%1e5;
% N_min = 1e4;%5e3;%1e5;
p_m = 1e-3;
% cat_c = 1; %cat_c_vec(p);% boolean variable to tell if there is catalysis of carbon fixation
% cat_p = 1;  % cat_p_vec(p);% boolean variable to tell if there is catalysis of polymerisation
%min_length_for_avg = 4;
%er = 0;%0.005;%[0,0.005,0.05,0.5];
er_c = er;%[0,0,0,0];%[0,0.005,0.05,0.5]; % error rate for copying
er_t = er;%[0,0.005,0.05,0.5]; % error rate for translation

% Initialise vectors
n_pops = cell(n,n_rec);
aa_pops = cell(n,n_rec);
k_avg_vec = zeros(n,n_gen);%zeros(1,n_gen);
avg_length_long_vec = zeros(n,n_gen);%zeros(1,n_gen);
N_m_vec_vec = zeros(n,n_gen);%zeros(1,n_gen);
p_c_vec_vec = zeros(n,n_gen);%zeros(1,n_gen);

% N_m_vec = cell(n,n_rec); 
% p_c_vec = zeros(length(p_p_min),t);
% total_mon_n_vec = zeros(length(p_p_min),t);
% avg_length_long_vec = zeros(length(p_p_min),t);
% num_pol_long_vec = zeros(length(p_p_min),t);

    
[t_mats_c] = error_fun(max_l,er_c); %Transition matrices (multiply these to the matrix that gives templates for copying/translating to represent error
[t_mats_t] = error_fun(max_l,er_t); %Transition matrices (multiply these to the matrix that gives templates for copying/translating to represent error

n_m_vec = zeros(1,n);
n_m_mat = zeros(n_gen,n);
div_per_gen = zeros(1,n_gen);


%% START SIMULATION
for g = 1:n_gen

    aa_pop_end = cell(1,n);
    n_pop_end = cell(1,n);
    k_avg = zeros(1,n);
    avg_length_long = zeros(1,n);
    N_m_vec = zeros(1,n);
    p_c_vec = zeros(1,n);
    
    for j = 1:n
        [n_pop_end{j}, aa_pop_end{j}, k_avg(j), total_mon_n, total_n_free_mon, total_aa_free_mon, total_n_mol,total_aa_mol,long_phobic, long_phylic, N_m_vec(j), p_c_vec(j), avg_length_long(j), num_pol_long] = model_fun(n_pop{j},aa_pop{j},alpha, beta, gamma, gamma2, epsilon,p_p_min,p_c_min,p_t_min,p_d,cat_c, cat_p,N_min, AA_min,p_m, min_length_for_avg,t_mats_c,t_mats_t,t,w);
    end
    
    
    %% Record information about the cells before they divide (only in the specified intervals)
    if rem(g,size_rec_int) == 0
        rec_n = g/size_rec_int;
        for h = 1:n
        n_pops{h,rec_n} = n_pop_end{h};
        aa_pops{h,rec_n} = aa_pop_end{h};
        end       
    end
    
    k_avg_vec(:,g) = k_avg;%k_avg_vec(g) = mean(k_avg);
    avg_length_long_vec(:,g) = avg_length_long;%avg_length_long_vec(g) = mean(avg_length_long);
    N_m_vec_vec(:,g) = N_m_vec;%N_m_vec_vec(g) = mean(N_m_vec);
    p_c_vec_vec(:,g) = p_c_vec;%p_c_vec_vec(g) = mean(p_c_vec);

    %% Division
    
    %Calculate number of nucleotide monomers per cell

    

    for k = 1:n
        n_m_vec(k) = total_mon_fun(n_pop_end{k},l_vec);
    end
    
    n_m_mat(g,:) = n_m_vec;
    disp(['gen',num2str(g),' ',num2str(find(n_m_vec >=threshold_size))]);

    aa_pop = aa_pop_end;
    n_pop = n_pop_end;
    div_per_gen(1,g) = sum(n_m_vec>=threshold_size);
    
if div_per_gen (1,g) >0
to_divide = find(n_m_vec>=threshold_size);

if length(to_divide)>= length(n_pop)/2
    disp('Warning! More than half the cells dividing at once!')
end

d1_aa = cell(1,length(to_divide));
d2_aa = cell(1,length(to_divide));
d1_n = cell(1,length(to_divide));
d2_n = cell(1,length(to_divide));
for i = 1:length(to_divide)

    [d1_aa{1,i},d2_aa{1,i},d1_n{1,i},d2_n{1,i}] = div_fun(aa_pop_end{1,to_divide(i)},n_pop_end{1,to_divide(i)},0.5,0.2);
    
end


% Replace mother cells with daughter 1
for j = 1:length(to_divide)
n_pop{to_divide(j)} = d1_n{j};
aa_pop{to_divide(j)} = d1_aa{j};
end


% Delete other random cell and replace with daughter 2

remaining_cells = 1:length(n_pop);
[~,b] = ismember(to_divide,remaining_cells);
remaining_cells(b) = [];
i_to_delete = randsample(length(remaining_cells),length(to_divide));
to_delete = remaining_cells(i_to_delete);

for k = 1:length(to_delete)
n_pop{to_delete(k)} = d2_n{k};
aa_pop{to_delete(k)} = d2_aa{k};
end

end



end

%n_gen = 1000
 %plot(1:n_gen, n_m_mat);
% 
% 
% x = movsum(div_per_gen,50, "Endpoints","discard");
% 
% plot(x);
% 
% 
