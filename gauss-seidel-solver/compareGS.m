% Artem Poddierohin
%function results = compareGS(A, B, tol, maxIter)
% compareGS   Run and compare forward vs. backward Gauss–Seidel
%
% Inputs:
%   A, B, tol, maxIter   – as above
%
% Output:
%   results  – struct with fields:
%     .X_gs, .iter_gs, .resHist_gs
%     .X_bgs, .iter_bgs, .resHist_bgs
%     .condA, .err_gs(relative), .err_bgs(relative)
%
% Structure:
% Part 1: Pre-computations
%   • condA = cond(A)  
%   • X_exact = A\B  
%
% Part 2: Run solvers
%   [X_gs,  iter_gs,  resHist_gs]  = gaussSeidel(A,B,tol,maxIter);
%   [X_bgs, iter_bgs, resHist_bgs] = backwardGaussSeidel(A,B,tol,maxIter);
%
% Part 3: 
%   • err_gs  = norm(X_gs  – X_exact)/norm(X_exact)
%   • err_bgs = norm(X_bgs – X_exact)/norm(X_exact)
%   • Package all into results struct
%end

function results = compareGS(A, B, tol, maxIter)

condA   = cond(A);        % Condition number
X_exact = A\B;            % MATLAB's reference solution

[X_gs,  iter_gs,  resHist_gs]  = gaussSeidel(A, B, tol, maxIter);
[X_bgs, iter_bgs, resHist_bgs] = backwardGaussSeidel(A, B, tol, maxIter);

err_gs  = norm(X_gs  - X_exact, 'fro') / norm(X_exact, 'fro');
err_bgs = norm(X_bgs - X_exact, 'fro') / norm(X_exact, 'fro');

results = struct('X_gs', X_gs, ...
                 'iter_gs', iter_gs, ...
                 'resHist_gs', resHist_gs, ...
                 'X_bgs', X_bgs, ...
                 'iter_bgs', iter_bgs, ...
                 'resHist_bgs', resHist_bgs, ...
                 'condA', condA, ...
                 'err_gs', err_gs, ...
                 'err_bgs', err_bgs);
end
