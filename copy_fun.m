function [new_copies] = copy_fun(N, N_e, n_philics,n_phobics)
%%% Copying function that takes the population of nucleotides (N), the
%%% population of templates, given an error rate, N_e, and matrices with
%%% the number of purines and pyrimidines per position, and returns a
%%% matrix with the new copies made.


% Check if there are enough monomers for everything
phobics_needed = N_e.*n_philics; %using the matrix for phobics here because of base pairing
philics_needed = N_e.*n_phobics; %using the matrix for phobics here because of base pairing

if sum(phobics_needed(2:end,:),'all') <= N(1,2) && sum(philics_needed(2:end,:),'all') <= N(1,1)
% if there are enough monomers for everything - proceed with copying
[pos_vec,~] = matrix2fvector(N_e,[], 'c'); %Make vector of positions to invert them because of base-pairing

[m,n] = size(N_e);
new_copies = zeros(m,n);
for o = 1:length(pos_vec)
    new_copies(pos_vec(o)) = new_copies(pos_vec(o))+1; %Turn vector above into a matrix that can be added to the population
end


else % else check which polymers can be made
%% Check what copies can be made given the amount of free purines in N

[pos_vec_pho,quant_vec_pho] = matrix2fvector(N_e,n_philics,'c');% Make frequency distribution from elements in the population to be copied, 
%note that the function already excludes monomers from being in the freq.
%dist. also note that the opposite matrix (philics vs phobics) is used to
%find the number of monomers because of base pairing

shuffled_i = randperm(length(quant_vec_pho));% Shuffle the vec of indices
cum_sum_pho = cumsum(quant_vec_pho(shuffled_i));% Make vector of cumulative sums, given shuffled indices
last_poss_pho = max(find(cum_sum_pho<=N(1,2)));% find last index that can be copied
poss_i_pho = shuffled_i(1:last_poss_pho);
pos_vec_new_pho = pos_vec_pho(poss_i_pho); % These are the positions that can be made given the number of purines/hydrophobics 


%% Check what copies can be made given the amount of pyrimidines in N

[pos_vec_phi,quant_vec_phi] = matrix2fvector(N_e,n_phobics,'c');% Make frequency distribution from elements in the population to be copied, 
%note that the function already excludes monomers from being 'in the freq. dist. also note that the opposite matrix (philics vs phobics) is used to
%find the number of monomers because of base pairing

cum_sum_phi = cumsum(quant_vec_phi(shuffled_i));% Make vector of cumulative sums, given shuffled indices
last_poss_phi = max(find(cum_sum_phi<=N(1,1)));% find last index that can be copied
poss_i_phi = shuffled_i(1:last_poss_phi);
pos_vec_new_phi = pos_vec_phi(poss_i_phi); % These are the positions that can be made given the number of purines/hydrophobics 

%% Find the copies for which there are enough monomers overall

pos_vec = intersectWithRepetitions(pos_vec_new_pho, pos_vec_new_phi);

%% Update the matrix
[m,n] = size(N_e);
new_copies = zeros(m,n);

for o = 1:length(pos_vec(:,1))
    new_copies(pos_vec(o,1)) = new_copies(pos_vec(o,1))+pos_vec(o,2);
end

    
end    