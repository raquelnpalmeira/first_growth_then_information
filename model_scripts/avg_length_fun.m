function [avg_length, num_pol] = avg_length_fun(pop,limit_length)

%Limit_length gives the max value of lengths to be excluded 

lengths = zeros(length(pop(:,1)), length(pop(1,:)));

for k = 1:length(pop(1,:))
for l = 1:length(pop(:,1))
    
    if k > l+1
        
        lengths(l,k) = 0;
    else
        lengths(l,k) = l;
    end
end
end




% Cut both matrices to exclude polymers with length smaller than the limit length 


% Multiply matrices element wise 

length_times_n = pop(limit_length+1:end,:).*lengths((limit_length+1):end,:);

avg_length = sum(length_times_n,"all") / sum(pop(limit_length+1:end,:),"all");

num_pol = sum(pop(limit_length+1:end,:),"all");

end