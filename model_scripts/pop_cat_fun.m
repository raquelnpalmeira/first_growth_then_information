function [k_avg,k_total] = pop_cat_fun(k_hs,gamma,epsilon,pop,hp_score,lengths,p_or_c)


%% Remove padding and monomers from matrices
%Initialising vector for indices
index = zeros(1,length(pop(:,1)));% do length -1 if excluding monomers

%Looping though lengths to see if there are any polymers of that length
for u = 1:length(pop(:,1))-1% if excluding monomers add %-1
    
    index(u+1) = sum(pop(u+1,:)>0)>0;%do u+1 if excluding monomers
    
end

%Making vector logical so we can index
index = logical(index);

%Indexing only lengths that have a polymer of that length
pop = pop(1:find(index,1, 'last'),:);
pop = pop(:,1:length(pop(:,1))+2); %this is because the number of columns should always be the number of rows +1, and in this case I have excluded mnomers

%Indexing the hp_score and lengths to match the population of peptides
hp_score = hp_score(1:find(index,1, 'last'),:);
hp_score = hp_score(:,1:length(hp_score(:,1))+2);
lengths = lengths(1:find(index,1, 'last'),:);
lengths = lengths(:,1:length(lengths(:,1))+2);



%% Loop through to calculate k values and sum

if sum(pop,'all') == 0 

    k_avg = 0;
    k_total = 0;
elseif sum(pop, 'all') > 0
    
    
% From hydrophobicity
kh = zeros(1,size(pop,1)*size(pop,2));
kl = zeros(1,size(pop,1)*size(pop,2));
k_total = 0;

if p_or_c == 1 % if 1, catalysis of polymerisation - better if hydrophilic (lower hp values)
    
    for i = 1:length(kh)
        %kh(i) = exp(hp_score(i)*k_hs) - exp(-k_hs);
        kh(i) = exp(-((hp_score(i)+k_hs)^2)/0.02)*0.0546; %(exp((hp_score*k_hs)))/1000;
        kl(i) = (1/(1+exp(-lengths(i)+gamma)))*epsilon;%kl(i) = kl_max * lengths(i)/(lengths(i)+ks);
        
        k_total = k_total + ((kh(i)*kl(i))*pop(i));
    end
    
else
    if p_or_c == 0 % if 0, catalysis of carbon fixation - better if hydrophobic (higher hp values)
        
        for i = 1:length(kh)
            %kh(i) = exp(hp_score(i)*(-k_hs)) - exp(-k_hs);
            kh(i) = exp(-((hp_score(i)+k_hs)^2)/0.02)*0.0546;%kh = (exp((-hp_score*k_hs)))/1000;
            kl(i) = (1/(1+exp(-lengths(i)+gamma)))*epsilon;%kl(i) = kl_max * lengths(i)/(lengths(i)+ks);
            
            k_total = k_total + ((kh(i)*kl(i))*pop(i));
        end
        
    else
        error('Invalid catalysis type')       
    end
end

k_avg = k_total/(sum(pop,'all'));
else
    error('Negative number of peptides');

end
end

