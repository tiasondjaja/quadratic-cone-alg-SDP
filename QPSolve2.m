function [ xOpt, val, solutionExists ] = QPSolve2( A, b, c, r, e )
% PROBLEM TO BE SOLVED:
% Compute the solution to QP(E, r):
% min trace(C*X)
% s.t. trace(A_i, X) = b_i
%      trace(E^{-1/2}*X*E^{-1/2})^2 - r^2*||E^{-1/2}*X*E^{-1/2}||^2 >= 0

%% PRELUDE
[ m, N ] = size( A ); %N = n*(n+1)/2
n = (-1+sqrt(1+8*N))/2;

%% Scale by E:
E = symvecinv(e);
sqrtE = sqrtm(E);

Atilde = zeros(m, N);
for i=1:m
    Ai = sqrtE*symvecinv(A(i, :)')*sqrtE;
    Atilde(i, :) = symvec(Ai)';
end
ctilde = symvec( sqrtE*symvecinv(c)*sqrtE );

idv = symvec(eye(n));

%% FRITZ JOHN

% 1. FRITZ JOHN optimality conditions to be solved -- linear part
Mat1 = -1/r^2 * ( (Atilde*idv)*(Atilde*idv)' / (n-r^2) - Atilde*Atilde' );
Mat2 = 1/r^2 * Atilde * ( idv*(idv'*ctilde) / (n-r^2) - ctilde );
Mat = [Mat1, Mat2];
rhs = b;

%  Find one possible Fritz-John point (one of many on a 1-dim'l subspace)
%    Use LDL' factorization
mu = 0;
[L, D] = ldl(Mat1);
y =  L' \ ( D\ (L\rhs) );
vmu = 1;
vy = -Mat1\Mat2;

res = (Atilde'*y - mu*ctilde);
vres = (Atilde'*vy - vmu*ctilde);

x = -1/r^2 * ( idv*(idv'*res) / (n-r^2) - res );
v = -1/r^2 * ( idv*(idv'*vres) / (n-r^2) - vres );

% 2. FRITZ JOHN optimality condition to be solved -- quadratic part
alpha = (v'*idv)^2 - r^2*(v'*v);
beta = 2*( (x'*idv)*(v'*idv) - r^2 * (x'*v) );
gamma = (x'*idv)^2 - r^2*(x'*x);
discriminant = beta^2 - 4*alpha*gamma;

solutionExists = true;
if ( discriminant < 0 )
    solutionExists = false;
else
    root1 = ( -beta - sqrt(discriminant) ) / ( 2*alpha );
    root2 = ( -beta + sqrt(discriminant) ) / ( 2*alpha );

    x1 = x + root1*v;
    x2 = x + root2*v;

    if ( (x1'*idv) >= 0 && (x2'*idv) >= 0 )
        
        if ctilde'*x1 <= ctilde'*x2
            xOpt = x1;
        else
            xOpt = x2;
        end
        val = ctilde'*xOpt;
        
    elseif ( (x1'*idv) >= 0 && (x2'*idv) < 0 )
        
        if ctilde'*x1 >= ctilde'*x2
            xOpt = x1;
            val = ctilde'*xOpt;
        else
            solutionExists = false;
            xOpt = zeros(N, 1);
        end
        
    elseif ( (x1'*idv) < 0 && (x2'*idv) >= 0 )
        
        if ctilde'*x1 <= ctilde'*x2
            xOpt = x2;
            val = ctilde'*xOpt;
        else
            solutionExists = false;
            xOpt = zeros(N, 1);
        end
        
    else
        solutionExists = false;
        xOpt = zeros(N, 1);
    end
end

if ~solutionExists
    val = inf;
    xOpt = x;
end

%% Unscale by inv(E):

XOpt = sqrtE*(symvecinv(xOpt)*sqrtE);
xOpt = symvec(XOpt);

end