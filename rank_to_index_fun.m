function [index1,index2] = rank_to_index_fun(rank,matrix) %% This exists as a matlab inbuilt function called index2sub, should try the code replacing this to avoid bugs 
size_rank = size(rank);
n_col_rank = size_rank(1);



if rank<= numel(matrix)
    size_matrix = size(matrix);
    index1 = ceil(rank/size_matrix(2));
    
    remainder = rem(rank,size_matrix(2));
    
    if remainder == 0
        
        index2 = size_matrix(2);
        
    else
        index2 = remainder;
        
    end
else
    error('Rank out of bounds');
end
end