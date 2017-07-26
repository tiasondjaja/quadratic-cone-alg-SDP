# quadratic-cone-alg-SDP
A MATLAB implementation of a quadratic-cone relaxation-based algorithm for semidefinite programming (https://arxiv.org/abs/1410.6734).


Files:
1. `test.m`

   A test script.  Calls `generateProblem` to randomly generate a test SDP and calls `QCRBASDP`, the main function, to solve the generated SDP.

2. [ A, b, c, seed ] =  generateProblem( m, n, a )
   
   Generates an SDP instance with variable X that is an nxn symmetric matrix, with m constraints.  
   
   + Inputs: 
     1. m = number of constraints
     2. n = number of rows/columns of the decision variable X
     3. a = specifies the seed for the first pseudorandom number generator.
   + Outputs:
     1. data for the SDP instance: A in R^(m x N), b in R^m, c in R^N (Note: N = n*(n+1)/2)
     2. seed = a record of all seeds for the randomly generated entries
     
2. `QCRBASDP.m`
   
   A "Quadratic Cone Relaxation-Based Algorithm for SDP".
   + Solves: Min trace(C*X), st. trace(A_i*X) = b, X in S^(nxn)
   + Inputs:
      1. data: A_i in S^(nxn), b in R^m, C in S^(nxn);
      2. initial iterate: E0 in S^(nxn);
      3. maxIt = max num of iterations;
      4. dualityGapBound = terminating condition
   + Outputs:   
      1. XOpt = final opt solution iterate;
      2. E = final center direction iterate;
      3. xOptVec = sequence of all symvec(XOpt)'s
      4. EVec = sequence of all E's
      5. val = final objective value
      6. exitFlag = 0: Optimal solution found; -1: Stops because reaches max number of iterations; -2: within the max number of iterations but duality gap is negative; -3: Else
    + Calls QPSolve2 to solve subproblems

3. `QPSolve2.m`
   
   + Solves QP(E, r): min trace(C*X), s.t. trace(A_i, X) = b_i; trace(E^{-1/2}*X*E^{-1/2})^2 - r^2*||E^{-1/2}*X*E^{-1/2}||^2 >= 0
   + Inputs: 
     1. data: A in R^(m x N), b in R^m, c in R^N (Note: N = n*(n+1)/2)
     2. cone-width parameter: r in (0, 1)
     3. initial iterate: e = a strictly feasible solution to the main SDP
   + Outputs:
     1. xOpt = the optimal solution
     2. val = the value of the optimal solution
     3. solutionExists = true if the input problem to QPSolve2 has an optimal solution
 
4. `initialize.m`
 
    + Finds a strictly feasible solution to the main SDP, by finding the analytic center of the feasible region.  Uses CVX
    + Inputs: Data: A in R^(m x N), b in R^m (Note: N = n*(n+1)/2 )
    + Output: e = a strictly feasible solution to the main SDP

5.  `symvec.m`

    + Converts an nxn symmetric matrix into a vector in R^N, where N = n*(n+1)/2
    + Input: A = an nxn symmetric matrix
    + Output: v = a vector in R^N
    
6.  `symvecinv.m`

    + Converts a vector in R^N into  an nxn symmetric matrix, where N = n*(n+1)/2
    + Input: v = a vector in R^N
    + Output: A = an nxn symmetric matrix

