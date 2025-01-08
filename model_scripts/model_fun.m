function [N,P, k_avg, total_mon_n, total_n_free_mon, total_aa_free_mon, total_n_mol,total_aa_mol, long_phobic, long_phylic, N_m_vec, p_c_vec, avg_length_long, num_pol_long] = model_fun(N,P, alpha,beta, gamma, gamma2, epsilon,p_p_min,p_c_min,p_t_min,p_d,cat_c, cat_p, N_min, AA_min, p_m, min_length_for_avg,t_mats_c,t_mats_t, t_max,w)

k_avg = zeros(1,t_max);
N_m_vec = zeros(1,t_max);
p_c_vec = zeros(1,t_max);
l_vec = 1:length(N(:,1));
total_mon_n = zeros(1,t_max);
total_n_free_mon = zeros(1,t_max);
total_aa_free_mon = zeros(1,t_max);
total_n_mol = zeros(1,t_max);
total_aa_mol = zeros(1,t_max);
long_phobic = zeros(1,t_max);
long_phylic = zeros(1,t_max);
avg_length_long = zeros(1,t_max);
num_pol_long = zeros(1,t_max);

%copy_index_vec = {};

for t = 1:t_max
%%  A matrix with hydrophobicity scores
hp_score = zeros(length(P(:,1)), length(P(1,:))); %subtract 1 from first dimention if excluding monomers 

for k = 1:length(P(1,:))
for l = 1:length(P(:,1)) %use 2:length to exclude monomers
    
    hp_score(l,k) = ((k-1)/l) - ((l-(k-1))/l); % if excluding monomers use the following: proportion_phobic(l-1,k) ;
    
end
end

smaller_than_1 = hp_score <= 1 ;
hp_score = hp_score .* smaller_than_1;

%% A matrix with the length of each polymer
lengths = zeros(length(P(:,1)), length(P(1,:)));

for k = 1:length(P(1,:))
for l = 1:length(P(:,1))
    
    if k > l+1
        
        lengths(l,k) = 0;
    else
        lengths(l,k) = l;
    end
end
end

%% Matrices with the number of hydrophibic and hydrophilic monomers needed for each position of the population 

n_phobics = zeros(length(N(:,1)), length(N(1,:)));

for i = 1:length(N(:,1))
    for j = 1:i+1
    n_phobics(i,j) = j-1;    
    end
end

n_philics = lengths - n_phobics;

%% Adding monomers

if sum(P,'all')>0
    [k_c,k_c_total] =  pop_cat_fun(alpha,gamma,epsilon,P,hp_score,lengths,0); % 0 for catalysis of carbon fixation
elseif sum(P,'all')==0
    k_c = 0; k_c_total= 0;
else
    error('negative sum of peptide numbers!');
end

if cat_c == 1 % If assuming hydrophobic peptides catalyse carbon fixation
    p_m_cat = p_m + (p_m*k_c_total);
    if p_m_cat >1
        p_m_cat = 1;
    end
    N_m = binornd(N_min, p_m_cat);
    AA_m = binornd(AA_min, p_m_cat);
elseif cat_c == 0
    N_m = binornd(N_min,p_m);%N_min*p_m;
    AA_m = binornd(AA_min,p_m);%AA_min*p_m;
else 
    error("Invalid cat_c value, must be either 0 or 1")
end


% Sample number of pyrimidines
add_Y = round(0.5 * N_m); % !!! This might create a bias as the number of pyrimidines will always be rounded up if the total number of nucleotides to be added is odd !!!
N(1,1) = N(1,1) + add_Y; % If making it more stochastic use: N_polymers(1,1) = binornd(N_min,0.5);
% Purines are the rest
N(1,2) = N(1,2) + (N_m - add_Y);

% Sample number of hydrophilic amino acids
add_B = round(0.5 * AA_m);
P(1,1) = P(1,1) + add_B; % If making it more stochastic use: P(1,1) = binornd(AA_min,0.5);
% Hydrophobic aas are the rest
P(1,2) = P(1,2) + (AA_m - add_B);

if rem(t,100)==0
disp(strcat("t = ", num2str(t)," ",num2str(N_m),"nucleotides added and ",num2str(AA_m),"amino acids added"));
end

monomers = sum(N.*lengths,'all');

%% Polymerising %For now only nucleotides "freely" polymerise

%If assuming hydrophilic amino acids catalyse polymerisation
% if cat_p == 1
%     k_p =  pop_cat_fun(alpha,kl_max,ks,P,hp_score,lengths,1); % 1 for catalysis of carbon fixation
% elseif cat_p ==0
%     k_p = 0;
% else
%     error("Invalid cat_p value, should be 1 or 0")
% end
% 
% p_p = p_p_min +(p_p_min*k_p);
% 
% if p_p >1 
%     p_p = 1;
% end

p_p = p_p_min;

% Sample pairs of molecules that include at least one monomer 
N_coords = pairs_to_pol_fun(N,p_p);

% Update matrix according to binding pairs
N = polfun(N_coords,N);


%% Copying 

total_n_free_mon(t) = sum(N(1,:));
total_aa_free_mon(t) = sum(P(1,:));
total_n_mol(t) = sum(N, 'all');
total_aa_mol(t) = sum(N, 'all');
[avg_length_long(t),num_pol_long(t)] = avg_length_fun(N, min_length_for_avg);

% If assuming hydrophilic amino acids catalyse polymerisation, calculate
% catalytic rate

[~,k_p_total] = pop_cat_fun(beta,gamma2,epsilon,P,hp_score,lengths,1); % 1 for catalysis of polymerisation

% Apply catalytic rate
if cat_p == 1
    p_c = p_c_min +(p_c_min*k_p_total);
    if p_c >1
        p_c = 1;
    end
elseif cat_p ==0
    p_c = p_c_min;
else
    error("Invalid cat_p value, should be 1 or 0")
end

if p_c >0.99
    disp('Warning, p_c > 0.99, normal approximation might give incorrect results.')
end

% Sample which will be templates given the new probability of copying
if p_c == 1
    N_temps = N;
elseif p_c <= 1 && p_c >=0 %Use normal approximation when values get too big
        N_temps = zeros(size(N));
        bino_val = N(N*p_c<=120);
        bino_sample = binornd(bino_val,p_c);
        N_temps(N*p_c<=120) = bino_sample;
        norm_val = N(N*p_c>120);
        norm_sample = round(normrnd(norm_val*p_c,sqrt(norm_val*p_c*(1-p_c))));
        norm_sample_t0 = max(0,norm_sample); %truncate so never negative
        norm_sample_t = min(norm_val,norm_sample_t0); %truncate so never bigger than the number of polymers available
        N_temps(N*p_c>120) = norm_sample_t;
else
    error("Invalid probability of copying");
end

 if sum(N_temps(2:end,:),'all')> 0 %Only update population if there are polymers to copy

% Find copies that would be made given error rate by multiplying by the
% transition matrices 
N_temps_e = zeros(size(N));
for l = 2:length(N(:,1)) %from 2 so that nothing is changed about monomers 
    round_m = round(N_temps(l,1:l+1)*t_mats_c{l});
    if sum(round_m) == sum(N_temps(l,1:l+1))
     N_temps_e(l,1:l+1) = round(N_temps(l,1:l+1)*t_mats_c{l});
    %disp(sum(N_temps(l,1:l+1),'all'));
    else
    N_temps_e(l,1:l+1) = floor(N_temps(l,1:l+1)*t_mats_c{l});
    left_over = sum(N_temps(l,1:l+1)) - sum(floor(N_temps(l,1:l+1)*t_mats_c{l}));
    %disp([sum(N_temps(l,1:l+1),'all'),left_over]);
    %pos_sample = datasample(1:l+1,left_over);
    [~,i_sort_decs] = sort(N_temps(l,1:l+1)*t_mats_c{l} - floor(N_temps(l,1:l+1)*t_mats_c{l}),'descend'); %get just the decimals by subtracting the integer from the whole number (i.e. number - floor(number)), then put them in descending order 
    i_sort_decs = i_sort_decs(1:left_over);
    for pos = i_sort_decs%pos_sample
    N_temps_e(l,pos)= N_temps_e(l,pos) +1;
    end
    end
end

%From the error corrected template make copies (only those for which we
%have enough monomers)
new_copies = copy_fun(N, N_temps_e, n_philics,n_phobics);

N = N+new_copies; %Add new copies to the population of nucleotides
N(1,1) = N(1,1) - sum(new_copies.*n_philics, 'all'); %subtract number of pyrimidines used, base pairing already accounted for in copy_fun
N(1,2) = N(1,2) - sum(new_copies.*n_phobics, 'all'); %subtract number of purines used

if isempty(N(N<0))==0
    error('Error! Negative number of polymers/monomers!');
end

end

long_phobic(t) = N(12,1);
long_phylic(t) = N(12,13);
%copy_index_vec{t} = pol_to_copy_index;


%% Translating

% If assuming hydrophilic amino acids catalyse polymerisation

[k_p,k_p_total] = pop_cat_fun(beta,gamma2,epsilon,P,hp_score,lengths,1); % 1 for catalysis of polymerisation


if cat_p == 1
    p_t = p_t_min +(p_t_min*k_p_total);
    if p_t >1
        p_t = 1;
    end
elseif cat_p ==0
    p_t = p_t_min;
else
    error("Invalid cat_p value, should be 1 or 0")
end

if p_t >0.99
    disp('Warning, p_t > 0.99, normal approximation might give incorrect results.')
end

% Sample which will be templates given the new probability of copying
% if p_t == 1
%     N_temps = N;
% elseif p_t < 1 && p_t >=0 % Use normal approximation for some values
%         N_temps = zeros(size(N));
%         bino_val = N(N*p_t<=120);
%         N_temps(N*p_t<=120) = binornd(bino_val,p_t);
%         norm_val = N(N*p_t>120);
%         N_temps(N*p_t>120) = round(normrnd(norm_val*p_t,norm_val*p_t*(1-p_t)));
% else
%     error("Invalid probability of translating");
% end

if p_t == 1
    N_temps = N;
elseif p_t <= 1 && p_t >=0 %Use normal approximation when values get too big
        N_temps = zeros(size(N));
        t_bino_val = N(N*p_t<=120);
        t_bino_sample = binornd(t_bino_val,p_t);
        N_temps(N*p_t<=120) = t_bino_sample;
        t_norm_val = N(N*p_t>120);
        t_norm_sample = round(normrnd(t_norm_val*p_t,sqrt(t_norm_val*p_t*(1-p_t))));
        t_norm_sample_t0 = max(0,t_norm_sample); %truncate so never negative
        t_norm_sample_t = min(t_norm_val,t_norm_sample_t0); %truncate so never bigger than the number of polymers available
        N_temps(N*p_t>120) = t_norm_sample_t;
elseif p_t == 0 
    N_temps = [];
else
    error("Invalid probability of translation");
end


if sum(N_temps(2:end,:),'all')> 0 %Only update population if there are polymers to copy

% Find copies that would be made given error rate by multiplying by the
% transition matrices 

% N_temps_e = zeros(size(N));
% for l = 2:length(N(:,1)) %from 2 so that nothing is changed about monomers
%     round_m = round(N_temps(l,1:l+1)*t_mats_t{l});
%     if sum(round_m) == sum(N_temps(l,1:l+1))
%         N_temps_e(l,1:l+1) = round(N_temps(l,1:l+1)*t_mats_t{l});
%     else
%         N_temps_e(l,1:l+1) = floor(N_temps(l,1:l+1)*t_mats_t{l});
%         left_over = sum(N_temps(l,1:l+1)) - sum(floor(N_temps(l,1:l+1)*t_mats_t{l}));
%         pos_sample = datasample(1:l+1,left_over);
%         for pos = pos_sample
%             N_temps_e(l,pos)= N_temps_e(l,pos) +1;
%         end
%     end
% end

N_temps_e = zeros(size(N));
for l = 2:length(N(:,1)) %from 2 so that nothing is changed about monomers 
    round_m = round(N_temps(l,1:l+1)*t_mats_t{l});
    if sum(round_m) == sum(N_temps(l,1:l+1))
     N_temps_e(l,1:l+1) = round(N_temps(l,1:l+1)*t_mats_t{l});
    %disp(sum(N_temps(l,1:l+1),'all'));
    else
    N_temps_e(l,1:l+1) = floor(N_temps(l,1:l+1)*t_mats_t{l});
    left_over = sum(N_temps(l,1:l+1)) - sum(floor(N_temps(l,1:l+1)*t_mats_t{l}));
    %disp([sum(N_temps(l,1:l+1),'all'),left_over]);
    %pos_sample = datasample(1:l+1,left_over);
    [~,i_sort_decs] = sort(N_temps(l,1:l+1)*t_mats_t{l} - floor(N_temps(l,1:l+1)*t_mats_t{l}),'descend'); %get just the decimals by subtracting the integer from the whole number (i.e. number - floor(number)), then put them in descending order 
    i_sort_decs = i_sort_decs(1:left_over);
    for pos = i_sort_decs%pos_sample
    N_temps_e(l,pos)= N_temps_e(l,pos) +1;
    end
    end
end

%From the error corrected template make copies (only those for which we
%have enough monomers)
new_transcripts = translation_fun(P, N_temps_e, n_philics,n_phobics); %note that templates are nucleotides, but the population P is given for monomers to be checked

P = P+new_transcripts; %Add new copies to the population of nucleotides

P(1,1) = P(1,1) - sum(new_transcripts.*n_philics, 'all'); %subtract number of pyrimidines used
P(1,2) = P(1,2) - sum(new_transcripts.*n_phobics, 'all'); %subtract number of pyrimidines used


end

%% Degrading polymers

%w = 1; %boolean variable to tell if decay is weighted on length or not (0 = not weighted, 1 = weighted)


if sum(N(2:end,:),'all')>0 %only degrade if there are polymers
% Sample which polymers to loose a monomer 
%N_decay_index = pol_to_copy_fun(N_polymers,); %maybe add option to not randomly permutate elements so to save computation time
N_decay_index = sample_pol_fun(N,p_d,lengths,w);
if isempty(N_decay_index) == 0 %Only update matrix if a polymer has been chosen
% Polymers decay (i.e. loose one monomer at a time)
[N, ~] = pol_decay_fun(N,N_decay_index); 
end
end 

if sum(P(2:end,:),'all')>0 %only degrade if there are polymers
% Sample which polymers to loose a monomer 
%AA_decay_index = pol_to_copy_fun(P,p_d); %maybe add option to not randomly permutate elements so to save computation time
AA_decay_index = sample_pol_fun(P,p_d,lengths,w);
if isempty(AA_decay_index) == 0 %Only update matrix if a polymer has been chosen
% Polymers decay (i.e. loose one monomer at a time)
[P, ~] = pol_decay_fun(P,AA_decay_index); 
end
end

% if rem(t,100)==0
% disp(strcat('Time ',num2str(t)));
% end

N_m_vec(t) = N_m;
p_c_vec(t) = p_c;
k_avg(t) = (k_p+k_c)/2;
total_mon_n(t) = total_mon_fun(N,l_vec);
avg_length_long(t) = avg_length_fun(N,min_length_for_avg);



%% The following commented section makes makes a time series of length and hydrophobicity 
% %% Average length over time
% 
% for j = 1:length(P(:,1))
%     
%     AA_n_per_length(i,j) = sum(P(j,1:j+1));
%     
% end
% AA_n_times_length = AA_n_per_length(i,:).*(1:length(AA_n_per_length(i,:)));
% AA_avg_length(h,i) = sum(AA_n_times_length)/sum(AA_n_per_length(i,:));
% 
% 
% for f = 1:length(N_polymers(:,1))
%     
%     N_n_per_length(i,f) = sum(N_polymers(f,1:f+1));
%     
% end
% N_n_times_length = N_n_per_length(i,:).*(1:length(N_n_per_length(i,:)));
% N_avg_length(h,i) = sum(N_n_times_length)/sum(N_n_per_length(i,:));
% 
% %% Average proportion of hydrophobic amino acids
% 
% AA_index_polymers = zeros(1,length(P(:,1))); % do length -1 if excluding monomers
% 
% for u = 1:length(P(:,1))% if excluding monomers add % -1
%     
%     AA_index_polymers(u) = sum(P(u,:)>0)>0; %do u+1 if excluding monomers
%     
% end
% AA_a = P .* hp_score; %use this to calculate over all nucleotide "strings"
% %AA_a = P(2:end,:) .* proportion_phobic; %use this to exclude
% %monomers
% AA_a = AA_a(logical(AA_index_polymers),:);
% AA_size_a = size(AA_a);
% AA_a = AA_a(:,1:AA_size_a(1)+1);
% 
% AA_avg_prop_phobic(h,i) = sum(AA_a,'all')/sum(P(1:end,:),'all'); %use 2:end if to exclude monomers
% 
% %% Average proportion of hydrophobic nucleotides
% 
% N_index_polymers = zeros(1,length(N_polymers(:,1)));% do length -1 if excluding monomers
% 
% for u = 1:length(N_polymers(:,1))% if excluding monomers add %-1
%     
%     N_index_polymers(u) = sum(N_polymers(u,:)>0)>0;%do u+1 if excluding monomers
%     
% end
% N_n = N_polymers .* hp_score; %use this to calculate over all nucleotide "strings"
% %N_a = N_polymers(2:end,:) .* proportion_phobic; %use this to exclude
% %monomers
% N_n = N_n(logical(N_index_polymers),:);
% N_size_n = size(N_n);
% N_n = N_n(:,1:N_size_n(1)+1);
% 
% N_avg_prop_phobic(h,i) = sum(N_n,'all')/sum(N_polymers(1:end,:),'all');  %use 2:end if to exclude monomers

%% Plot histogram over time 
% freq_data = repelem(1:length(sum(N_polymers,2)),sum(N_polymers,2)); %fequency data for histogram
% figure(1)
% if rem(i,t/4) == 0 
%     subplot(2,2,(i/(t/4)));
%     histogram(freq_data,'FaceColor',colour{h},'EdgeColor',colour{h},'Normalization','probability');%using histogram because it supports transparency in the bars - https://uk.mathworks.com/matlabcentral/answers/213991-matlab-put-transparancy-on-a-bar-plot
%     title(['t=',num2str(i)]);
%     xlabel('length')
%     ylabel('frequency')
%     set(gca,'FontSize',20)
% end

 %% Plot an area 
% if rem(i,t/4) == 0 
%     j = j+1;
%     subplot(2,2,(i/(t/4)));
%     a = area(sum(N_polymers,2)/sum(N_polymers,"all"),'FaceColor',colour{h},'EdgeColor',colour{h});
%     a.Facegamma = 0.2;
%     title(['t=',num2str(i)]);
%     xlim([1 length(N_polymers(:,1))])
%     xlabel('length')
%     ylabel('frequency')
%     set(gca,'FontSize',18)
% end

%% Plot a bar chart 
% 
% if rem(i,t/4) == 0 
%     subplot(2,2,(i/(t/4)));
%     bar(sum(N_polymers,2)/sum(N_polymers,"all"),weight(h),'FaceColor',colour{h});
%     title(['t=',num2str(i)]);
%     xlabel('length')
%     ylabel('frequency')
%     set(gca,'FontSize',20)
% end
end
end
