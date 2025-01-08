function [index] = sample_pol_fun(matrix,p,lengths,w)

matrix_no_monomers = matrix;
matrix_no_monomers(1,:) = [];
matrix_no_monomers_trans = matrix_no_monomers';
matrix_no_monomers_vec = matrix_no_monomers_trans(:)'; %Make a vector with all the elements in matrix
matrix_no_monomers_dist = repelem(1:numel(matrix_no_monomers_vec),matrix_no_monomers_vec); %Make distribution to sample from

if w == 0 % Without weights based on length
    matrix2copy = datasample(matrix_no_monomers_dist, round(p*length(matrix_no_monomers_dist)),'Replace',false);
elseif w ==1 % With weights based on length
    lengths_no_mon = lengths;
    lengths_no_mon(1,:) = [];
    lengths_no_mon_trans = lengths_no_mon';
    if length(matrix_no_monomers_dist)>1
        weights = lengths_no_mon_trans(matrix_no_monomers_dist); %weights are the corresponding lengths
        matrix2copy = datasample(matrix_no_monomers_dist,round(p*length(matrix_no_monomers_dist)),'Replace',false,'Weights',weights);
    elseif length(matrix_no_monomers_dist)==1
        coin_toss = binornd(1,p);
        if coin_toss == 1
            matrix2copy = matrix_no_monomers_dist;
        else
            matrix2copy = [];
        end
    else
        error('no polymers to sample from!')
    end
else
    error('invalid weighting option')
end

if isempty(matrix2copy) == 1
    index = [];
else
    index = zeros(2,length(matrix2copy));
    [index(1,:), index(2,:)] = rank_to_index_fun(matrix2copy,matrix_no_monomers); % can replace this by inbuild function sub2ind
    index = index';
    index(:,1) = index(:,1)+1;
    %ranks = sub2ind(size(matrix), index(:,1)', index(:,2)');
end
end

