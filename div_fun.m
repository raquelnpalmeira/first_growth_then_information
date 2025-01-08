%% Division function
function[daughter1_aa,daughter2_aa,daughter1_n,daughter2_n] = div_fun(mother_aa,mother_n,mu,sigma)

    pd = makedist('Normal','mu',mu,'sigma',sigma); %make distribution to sample from
    t_pd = truncate(pd,0,1); %truncate it between 1 and 0


 %% FOR AMINO ACIDS
        % Find ouy how many types of polymers exist in the cell
        I_pol_aa = find(mother_aa ~= 0); % find the indexes for the positions that have polymers

        % find out how many types of polymers there are in the cell
        aa_types_pol = length(I_pol_aa);


        % sample proportions that will go into each daughter cell
        for_cell_1_aa = random(t_pd, aa_types_pol,1);
        for_cell_2_aa = 1 - for_cell_1_aa;


        % assign polymers to each daughter cell
        daughter1_aa = mother_aa;
        daughter2_aa = mother_aa;
        a  = round(mother_aa(I_pol_aa).*for_cell_1_aa);
        b = round(mother_aa(I_pol_aa).*for_cell_2_aa);

        if sum(a+b ~= mother_aa(I_pol_aa))>0
            error('Sum of daughter cell molecules more than molecules in mother cell')
        end

        daughter1_aa(I_pol_aa) = a;
        daughter2_aa(I_pol_aa) = b;

        if sum(daughter1_aa+daughter2_aa ~= mother_aa, 'all')
            error('Sum of daughter cell molecules more than molecules in mother cell')
        end

  %% FOR NUCLEOTIDES
        % Find ouy how many types of polymers exist in the cell

        I_pol_n = find(mother_n ~= 0); % find the indexes for the positions that have polymers

        % find out how many types of polymers there are in the cell
        n_types_pol = length(I_pol_n);


        % sample proportions that will go into each daughter cell
        for_cell_1_n = random(t_pd, n_types_pol,1);
        for_cell_2_n = 1 - for_cell_1_n;


        % assign polymers to each daughter cell
        daughter1_n = mother_n;
        daughter2_n = mother_n;
        c  = round(mother_n(I_pol_n).*for_cell_1_n);
        d = round(mother_n(I_pol_n).*for_cell_2_n);

        if sum(c+d ~= mother_n(I_pol_n))>0
            error('Sum of daughter cell nucleotides more than molecules in mother cell')
        end

        daughter1_n(I_pol_n) = c;
        daughter2_n(I_pol_n) = d;

        if sum(daughter1_n+daughter2_n ~= mother_n, 'all')
            error('Sum of daughter cell necleotides more than molecules in mother cell')
        end

    end