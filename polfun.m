function matrix = polfun(pairs,matrix)

[nrow_pairs,~] = size(pairs);

for i = 1:nrow_pairs
    x = pairs(i,1);
    y = pairs(i,2);
    a = pairs(i,3);
    b = pairs(i,4);
   
    matrix(a+x, b+(y-1)) = matrix(a+x, b+(y-1))+1;
    matrix(a,b) = matrix(a,b)-1;
    matrix(x,y) = matrix(x,y)-1;
    
end

end