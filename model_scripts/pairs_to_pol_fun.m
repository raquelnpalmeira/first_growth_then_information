function [coords] = pairs_to_pol_fun(matrix,p)
%Sample pairs of molecules 
matrix_trans = matrix';
matrix_vec = matrix_trans(:)'; %Make a vector with all the elements in matrix
matrix_dist = repelem(1:numel(matrix_vec),matrix_vec); %Make distrubution to sample from
matrix_perm = matrix_dist(randperm(length(matrix_dist))); %Randomly permutate the distribution

%Make random pairs
%i=1;
pairs = zeros(floor(length(matrix_perm)/2), 2);
% while length(matrix_perm) >= 2
%     pairs(i,:) = matrix_perm(1:2);
%     matrix_perm(1:2) = [];
%     i = i+1;
% end

for i = 1:floor(length(matrix_perm)/2)
    pairs(i,1) = matrix_perm((2*i)-1);
    pairs(i,2) = matrix_perm(2*i);
end



%pairs(:,1) = matrix_perm(1:2:end);
%pairs(:,2) = matrix_perm(2:2:end);

%Select only pairs that have at least one monomer (i.e. because we're
%considering polymerisation only happens with one monomer being added at a
%time
pairs_inc_monomers = pairs(or(any(pairs == 1, 2), any(pairs == 2, 2)),:);
size_pairs_inc_mon = size(pairs_inc_monomers);
% pairs_inc_mon_suc = randsample(1:size_pairs_inc_mon(1), round(p*size_pairs_inc_mon(1))); %Only some pairs are succesful depending on the probability
% pairs_inc_monomers = pairs_inc_monomers(pairs_inc_mon_suc,:);
pairs_inc_mon_suc = rand(size_pairs_inc_mon(1),1);
pairs_inc_monomers = pairs_inc_monomers(pairs_inc_mon_suc<p,:);

%Find indexes for the selected pairs
size_pairs = size(pairs_inc_monomers);
coords = zeros(size_pairs(1),size_pairs(2)*2);

for j = 1:size_pairs(1)
    for i = 1:size_pairs(2)
    
        [index1,index2] = rank_to_index_fun(pairs_inc_monomers(j,i), matrix);
        coords(j,2*i-1:2*i) = [index1,index2];
    end 
end

end



%% maybe sample a monomer at a time and pair it with another molecule