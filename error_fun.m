% Create transition matrices
function [t_matrices] = error_fun(max_l,er)


t_matrices = cell(1,max_l);


for l = 1:max_l
%% Initialising matrix

t_matrix = zeros(l+1,l+1);


for i = 0:l
    
    for j = 0:l
        
        k_vec = 0:j; k_vec = k_vec(k_vec<=i); k_vec = k_vec(k_vec<=j); k_vec = k_vec(k_vec>=i+j-l);
        sum_probs = 0;
        
        for k = k_vec
            
            perms_right = nchoosek(i,k);
            prob_right = perms_right*((1-er)^k)*(er^(i-k));
            
            perms_wrong = nchoosek(l-i,j-k);
            prob_wrong = perms_wrong*(er^(j-k))*((1-er)^((l-i)-(j-k)));
            
            new_prob = prob_right * prob_wrong;
            sum_probs = sum_probs+new_prob;
        end
        t_matrix(i+1,j+1) = sum_probs;
    end
end

tol = eps("single"); %very small number as tolerance %% I addedd this because it was giving me 1 is not equal to 1

if sum(abs(sum(t_matrix,2)-1))>tol %comparing the difference to a very small number
error(['for er = ',num2str(er),' sum of probabilities in transition matrix different to 1'])
end

t_matrices{l} = t_matrix;
end