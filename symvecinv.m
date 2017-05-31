function [ A ] = symvecinv( v )
%SYMVECINV takes an N = n*(n+1)/2-column vector v and
%   outputs the corresponding nxn symmetric matrix A

    N = size(v);
    n = (sqrt(1+8*N)-1)/2; % quadratic formula to solve for n
    
    for i = 1:n
        for j = i:n
           if i ~= j
               A(i, j) = v(i + j*(j-1)/2)/sqrt(2);
               A(j, i) = A(i, j);
           else
               A(i, i) = v(i + i*(i-1)/2);
           end
        end 
    end

end

