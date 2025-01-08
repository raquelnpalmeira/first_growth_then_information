function [matrix, R_decay_index, Y_decay_index] = pol_decay_fun(matrix,decay_index)

% Sample from distribution made from polymer whether the monomer lost will
%be purine or pyrimidine
p = (decay_index(:,2) - 1)./decay_index(:,1); %proportion of the string made of purines/hydrophobic aas 
p = p(:,1);
random_vector = rand(size(p));
%n = repelem(1,length(p))';
R_boolean = p>random_vector;%binornd(n,p);
Y_boolean = p<random_vector;%abs(R_boolean-1);   

% Add monomers that have "come out"
matrix(1,1) = matrix(1,1) + sum(Y_boolean);
matrix(1,2) = matrix(1,2) + sum(R_boolean);


% Update polymers (if they were dimers in the first place, this will update
% monomer numbers too)

R_logical = logical(R_boolean);
Y_logical = logical(Y_boolean);

R_decay_index = decay_index(R_logical,:);
Y_decay_index = decay_index(Y_logical,:);



for i = 1:length(R_decay_index(:,1))
   
    x = R_decay_index(i,1);
    y = R_decay_index(i,2);
    
    matrix(x-1,y-1) = matrix(x-1,y-1) +1;
    matrix(x,y) = matrix(x,y) -1;

end

for i = 1:length(Y_decay_index(:,1))
   
    x = Y_decay_index(i,1);
    y = Y_decay_index(i,2);
    
    matrix(x-1,y) = matrix(x-1,y) +1;
    matrix(x,y) = matrix(x,y) -1;

end

end