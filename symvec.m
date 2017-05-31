function [ v ] = symvec( A )
%SYMVEC takes a symmetric nxn matrix A and 
%   outputs an n*(n+1)/2-column vector v containing
%   transformed columns of the upper triangular part of A
%   (transformed to preserve inner product)

[~, n] = size(A);
v = zeros(n*(n+1)/2, 1);

for i = 1:n
    
    % if j = i
    v( i + i*(i-1)/2 ) = A(i, i);
    
    % if j\neq i
    for j = i+1:n
        v( i + j*(j-1)/2 ) = sqrt(2)*A(i, j);
    end
end

end

