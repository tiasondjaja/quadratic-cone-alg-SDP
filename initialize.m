function [ e ] = initialize( A, b )

    [ m, N ] = size(A);
    n = (-1+sqrt(1+8*N))/2;
    
    cvx_begin sdp
        variable E(n, n) symmetric;
        minimize -log_det(E)
        subject to
            for i=1:m
                trace(symvecinv(A(i, :)') * E) == b(i);
            end
            E == semidefinite(n)
    cvx_end
    
    e = symvec(E);
end

