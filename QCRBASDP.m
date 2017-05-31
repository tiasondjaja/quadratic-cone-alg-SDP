function [ XOpt, E, xOptVec, eVec, val, exitFlag, it, dualityGap, tVec, solutionExists, longStepItTotal ] = ...
    QCRBASDP( A, b, c, r, e0, maxIt, dualityGapBound, longStep0 )
% QCRBASDPv1 - Quadratic Cone Relaxation-Based Algorithm for SDP
%             solves: Min trace(C*X), st. trace(A_i*X) = b, X in S^(nxn)
% Input:    data: A_i in S^(nxn), b in R^m, C in S^(nxn);
%           initial iterate: E0 in S^(nxn);
%           maxIt = max num of iterations;
%           dualityGapBound = terminating condition
% Output:   XOpt = final opt solution iterate;
%           E = final center direction iterate;
%           xOptVec = sequence of all symvec(XOpt)'s
%           EVec = sequence of all E's
%           val = final objective value
%           exitFlag =
%                0: Optimal solution found
%               -1: Stops because reaches max number of iterations
%               -2: within the max number of iterations
%                   but duality gap is negative
%               -3: Else

    %% Variable declaration-initialization
    dualityGap = 1000;
    [~, N] = size(A);
    
    E = symvecinv(e0);
    e = e0;
    val = c'*e;
    
    xOptVec = zeros(0, N);
    eVec = zeros(0, N);
    tVec = zeros(maxIt, 1);
    
    %% Checking that E is in Swath(r)
    
     [ xOutput, ~, solutionExists ] = QPSolve2( A, b, c, r, e );
     
     if ( ~solutionExists )
         disp( 'Initial vector is not in Swath(r).  Abort!' );
         pause;
         
     end

    %% Main Loop
    
    it=1;
    longStep = longStep0;
    longStepIt = 0;
    longStepItTotal = 0;
    while ( it < maxIt && abs(dualityGap) > dualityGapBound && solutionExists)
        
        
        % Iteration i
        xOpt = xOutput;
        xOptVec = [ xOptVec; xOpt' ];
        eVec = [ eVec; e' ];
        
        % Determining search direction for computing next iterate
        XOpt = symvecinv( xOpt );
        val = c'*e;
        d = xOpt - e; % search direction
        dualityGap = -c'*d; %duality gap at iteration i
        
        % Determining step size for computing next iterate
        E = symvecinv( e );
        sqrtE = sqrtm( E );
        
        % short step
        t1 = r/norm( sqrtE\(XOpt/sqrtE), 'fro' )/2;
        tshortstep = t1/( 1+t1 );
        
        if (~longStep) % use step size prescribed in paper
            t = tshortstep;
            
            [ xOutput, ~, solutionExists ] = QPSolve2( A, b, c, r, e + t*d );
        else
            % long-step
            D = symvecinv( d );
            Dtilde = sqrtE\( D/sqrtE );
            eigenvalues = eig( Dtilde );
            negEigenvalues = eigenvalues( eigenvalues < 0 );
            tlongstep = 1/( -max( negEigenvalues ) )/2;
            t = max( tlongstep, tshortstep );
            
            [ xOutput, ~, solutionExists ] = QPSolve2( A, b, c, r, e + t*d );
        end
        
        % Check if long step produces iterate in Swath(r)
        while ( ~solutionExists && longStepIt < 4 && longStep ) 
            longStepIt = longStepIt + 1;
            t = t/2; % If not, 
            [ xOutput, ~, solutionExists ] = QPSolve2( A, b, c, r, e + t*d );
            
            if ( longStepIt == 4 && ~solutionExists )
                t = tshortstep;
                [ xOutput, ~, solutionExists ] = QPSolve2( A, b, c, r, e + t*d );
            end
        end
        longStepItTotal = longStepItTotal + longStepIt;
        longStepIt = 0; % reset value of longStepIt to zero
        
        tVec( it ) = t; % record the step length
        
        % Iteration i + 1
        % New iterate
        it = it + 1;
        e = e + t*d; % new iterate
        
    end    
    
    % Exit flag
    if ( dualityGap <= dualityGapBound && dualityGap >= 0 && solutionExists )
        exitFlag = 0; % stops because duality gap is below bound
    elseif ( it >= maxIt && dualityGap >= 0 )
        exitFlag = -1; % stops because reaches max number of iterations
    elseif ( it <= maxIt && dualityGap < 0  )
        exitFlag = -2; % within the max number of iterations
                      % but duality gap is negative
    else %else (e.g. E not in Swath(r) )
        exitFlag = -3;
    end
    
end
