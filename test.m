% test.m
% Randomized test case for QCRBASDPv1
% - Quadratic Cone Relaxation-Based Algorithm for SDP, v.1.
% PROBLEM TO BE SOLVED:
% Minimize trace(C*X), s.t. trace(A_i*X) = b_i (for i=1, ..., m) , X psd

clear all; 
close all;

%% INPUT
m = 10; %number of constraints
n = 10; %number of rows/cols of matrix variable X
a = 5;  %seed for first entry
[ A, b, c, seed ] = generateProblem( m, n, a );
C = symvecinv(c);

%% Solve SDP using the QCRBASDP
% PARAMETERS
r = 0.9;       %cone-width parameter
maxIt = 150;   %maximum number of iterations
dualityGapBound = 1e-5;
longStep = false;  %long steps vs. short steps

%% Determine initial center direction e for QCRBASDP
e0 = initialize( A, b );
e = e0;

%% Run QCRBASDP
[ XOpt, EOpt, xOptVec, eOptVec, val, exitFlag, it, dualityGap, tVec, solutionExists, itTotal ] = ...
    QCRBASDP( A, b, c, r, e, maxIt, dualityGapBound, longStep );

%% Solve SDP using CVX (for comparison purposes)
cvx_begin sdp
    variable X(n, n)
    minimize trace(C*X)
    for i=1:m
        Ai = symvecinv(A(i, :)');
        trace(Ai*X) == b(i);
    end
    X == semidefinite(n)
cvx_end
Xcvx = X;


%% Display output on command window
disp( ['Objective Value: ', num2str(trace(C*EOpt)) ] );
disp( ['Objective Value of CVX output: ', num2str(trace(C*Xcvx)) ] );
disp( ['Relative error (norm, compared to CVX output): ', num2str(norm(EOpt-Xcvx, 'fro')/norm(Xcvx, 'fro'))] )
disp( ['dualityGap: ', num2str(dualityGap)] );
disp( ['Number of iterations: ', num2str(it)]);
disp( ['exitFlag: ', num2str(exitFlag) ] );


%% Write output to text file
printResult = false;  %if true, results will be saved in a text file

if (printResult)
    fileID = fopen('test-output.txt','a');
    fprintf(fileID, '---------------\r\n');
    fprintf(fileID, '---------------\r\n');
    fprintf(fileID, 'Date: %s \r\n', date);
    fprintf(fileID,'m = %d, n = %d \r\n', m, n);
    fprintf(fileID, 'Seed vector: ');
    fprintf(fileID, '%d ', seed);
    fprintf(fileID,' \r\n');
    fprintf(fileID, 'r = %f \r\n', r);
    fprintf(fileID, 'long step? %d \r\n', longStep);
    fprintf(fileID, 'maxIt: %d \r\n', maxIt);
    fprintf(fileID, 'duality gap bound: %5e \r\n', dualityGapBound);
    fprintf(fileID, '---------------\r\n');
    fprintf(fileID,'Objective Value: %5e \r\n', trace(C*EOpt) );
    fprintf(fileID,'Objective Value of CVX output: %5e \r\n', trace(C*Xcvx) );
    fprintf(fileID,'Relative error, compared to CVX output (frobenius norm): %5e \r\n', norm(EOpt-Xcvx, 'fro')/norm(Xcvx, 'fro'));
    fprintf(fileID,'Duality gap: %5e \r\n', dualityGap);
    fprintf(fileID,'Number of iterations: %d \r\n', it);
    fprintf(fileID,'exitFlag: %d \r\n\r\n', exitFlag);
    fclose(fileID);
end
