function common_and_counts = intersectWithRepetitions(A, B)
    % INTERSECTWITHREPETITIONS Find elements common to both vectors A and B, including repetitions.
    %   common = INTERSECTWITHREPETITIONS(A, B) returns a vector containing the elements that are common
    %   to both input vectors A and B, including any repetitions. The output is sorted in ascending order.
    %
    %   Example:
    %       A = [1, 2, 2, 3];
    %       B = [2, 3, 3, 4];
    %       common_and_counts = INTERSECTWITHREPETITIONS(A, B)
    %       common_and_counts =
    %           2     1
    %           3     1
    %

    % Find the elements common to both vectors
    common = intersect(A, B);

    % Calculate the number of times each common element appears in both vectors
    countsA = zeros(1,length(common));
    countsB = zeros(1,length(common));

    for i = 1:length(common)
    countsA(i) = sum(A == common(i));
    countsB(i) = sum(B == common(i));
    end

    min_count = zeros(1,length(common));
    for j = 1:length(common)
    min_count(j) = min(countsA(j),countsB(j));
    end

    % Sort the output in ascending order based on the elements
    common_and_counts = [common',min_count'];
end
