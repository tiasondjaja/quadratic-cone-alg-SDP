function [ A, b, c, seed ] =  generateProblem( m, n, a )
%GENERATEPROBLEM generates inputs to an SDP problem with 
%     m constraints, with
%     the cone of nxn symmetric positive semidefinite matrices
%     a is the seed for the first randomly generated entry.

    N = n*(n+1)/2;
    % set seed for generating input:
    seed = zeros(m+2, 1);
    seed(1) = a;
    rng(seed(1));
    E0 = tril(rand(n, n));
    E0 = E0*E0';
    e0 = symvec(E0);
    A = zeros(m, N);
    for i=1:m
        seed(i+1) = i*n^2 + a + 1;
        rng( seed(i+1) );
        L = tril(rand(n, n));
        A(i, :) = symvec(L*L')';
    end % initial derivative direction
    b = A*e0; % rhs vector of the linear constraint
    seed(m+2) = (m+1)*n^2 + a + 1;
    rng( seed(m+2) );
    C = tril(rand(n, n));
    C = C*C'; % objective vector
    c = symvec(C);

end

