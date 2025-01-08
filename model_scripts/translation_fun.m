function [new_transcripts] = translation_fun(A, A_e, n_philics,n_phobics)
%%% Copying function that takes the population of nucleotides (N), the
%%% population of templates, given an error rate, N_e, and matrices with
%%% the number of purines and pyrimidines per position, and returns a
%%% matrix with the new copies made.


% Check if there are enough monomers for everything
phobics_needed = A_e.*n_phobics; 
philics_needed = A_e.*n_philics; 

if sum(phobics_needed(2:end,:),'all') <= A(1,2) && sum(philics_needed(2:end,:),'all') <= A(1,1)
% if there are enough monomers for everything - proceed with copying
[pos_vec,~] = matrix2fvector(A_e,[], 't'); %Make vector of positions

[m,n] = size(A_e);
new_transcripts = zeros(m,n);
for o = 1:length(pos_vec)
    new_transcripts(pos_vec(o)) = new_transcripts(pos_vec(o))+1; %Turn vector above into a matrix that can be added to the population
end


else % else check which polymers can be made
%% Check what copies can be made given the amount of free purines in N

[pos_vec_pho,quant_vec_pho] = matrix2fvector(A_e,n_phobics,'t');% Make frequency distribution from elements in the population to be copied, 
%note that the function already excludes monomers from being in the freq.
%dist. also note that the opposite matrix (philics vs phobics) is used to
%find the number of monomers because of base pairing

shuffled_i = randperm(length(quant_vec_pho));% Shuffle the vec of indices
cum_sum_pho = cumsum(quant_vec_pho(shuffled_i));% Make vector of cumulative sums, given shuffled indices
last_poss_pho = min(find(cum_sum_pho>=A(1,2)))-1;% find last index that can be copied
poss_i_pho = shuffled_i(1:last_poss_pho);
pos_vec_new_pho = pos_vec_pho(poss_i_pho); % These are the positions that can be made given the number of purines/hydrophobics 


%% Check what copies can be made given the amount of pyrimidines in N

[pos_vec_phi,quant_vec_phi] = matrix2fvector(A_e,n_philics,'t');% Make frequency distribution from elements in the population to be copied, 
%note that the function already excludes monomers from being 'in the freq. dist. also note that the opposite matrix (philics vs phobics) is used to
%find the number of monomers because of base pairing

cum_sum_phi = cumsum(quant_vec_phi(shuffled_i));% Make vector of cumulative sums, given shuffled indices
last_poss_phi = min(find(cum_sum_phi>=A(1,1)))-1;% find last index that can be copied
poss_i_phi = shuffled_i(1:last_poss_phi);
pos_vec_new_phi = pos_vec_phi(poss_i_phi); % These are the positions that can be made given the number of purines/hydrophobics 

%% Find the copies for which there are enough monomers overall

pos_vec = intersectWithRepetitions(pos_vec_new_pho, pos_vec_new_phi);

%% Update the matrix
[m,n] = size(A_e);
new_transcripts = zeros(m,n);

for o = 1:length(pos_vec(:,1))
    new_transcripts(pos_vec(o,1)) = new_transcripts(pos_vec(o,1))+pos_vec(o,2);
end

    
end    