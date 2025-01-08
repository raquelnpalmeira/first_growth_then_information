function [pos_vec,quant_vec] = matrix2fvector(A,mon_mat, type)

% Vector with positions in the matrix
[m, n] = size(A); % Get the dimentions of the matrix
size_A = m*n;
pos_vec = []; % Initialize the output vector

for j = 1:size_A
    if j ~= 1 && j ~= m+1 % This excludes monomers, since we don't want them being copied
        %             num_repeats = A(j); % Get the number of times to repeat the position
        %             rep_pos = repelem(j,num_repeats); % Repeat position j for the number of elements in the matrix
        rep_pos = j*ones(1,A(j));
        pos_vec = [pos_vec, rep_pos]; % Add the repeated positions to the output vector

    end
end

if type == 'c'
    %Transform pos_vec so that it reflects base-pairing %%
    [r,c]= ind2sub(size(A), pos_vec);
    c = r-(c-1)+1;
    pos_vec = sub2ind(size(A),r,c);
elseif type == 't'
else
    error('Incorrect type')
end


if isempty(mon_mat) == 1
    quant_vec = zeros(1,sum(A,'all'));

else
    % Vector with the number of purines or pyrimidines needed
    num_repeats = 0;

    quant_vec = []; % Initialize the output vector
    for j = 1:size_A
        if j ~= 1 && j ~= m+1 % This excludes monomers, since we don't want them being copied

            %             num_repeats = A(j); % Get the number of times to repeat the position
            %             num_mon = mon_mat(j);
            %             rep_mon = repelem(num_mon,num_repeats); % Repeat position j for the number of elements in the matrix
            rep_mon = mon_mat(j)*ones(1,A(j));
            quant_vec = [quant_vec, rep_mon]; % Add the repeated positions to the output vector
        end

    end
end
end
