function [total_mon] = total_mon_fun(x, l_vec)

% add across rows
sum_rows = sum(x,2);

% create vector with lengths 
l_vec = 1:length(x(:,1));

% multiply sums by lengths
total_mon_rows =  sum_rows .* l_vec';

% find total amound of monomers in the system
total_mon = (sum(total_mon_rows));

end