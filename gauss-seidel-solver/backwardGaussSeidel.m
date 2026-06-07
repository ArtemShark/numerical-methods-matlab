% Artem Poddierohin
%function [X, iter, resHistory] = backwardGaussSeidel(A, B, tol, maxIter)
% backwardGaussSeidel   Solve A*X = B using reverse-order Gauss–Seidel
%
% Inputs:
%   A        – n×n coefficient matrix, must be nonsingular so that the method makes sense
%   B        – n×m right-hand side
%   tol      – convergence tolerance on relative change
%   maxIter  – maximum iterations allowed 
%
% Outputs:
%   X            – n×m solution matrix 
%   iter         – number of iterations performed before stopping
%   resHistory   – vector of ||B–A*X_k|| for k=0:iter
%
% Structure:
% Part 1: Initialization
%   • Validate dimensions, initialize X, 
%   compute the resHistory(1) = norm(B - A*X)
%   This tells how far your zero-guess is from satisfying the system.
%
% Part 2: Iterative loop
%   for k = 1:maxIter Each iteration k runs one backward sweep across all rows
%     for i = n:-1:1                 % backward sweep
%       X(i,:) = (B(i,:) 
%                – A(i,1:i-1)*X(1:i-1,:) uses the newly computed entries(since we’re sweeping downward, those correspond to indices j>i)
%                – A(i,i+1:n)*X(i+1:n,:)) / A(i,i); uses the old entries from before this sweep (for j < i)
%     end
%     resHistory(k+1) = norm(B–A*X);  
%     if norm(X–X_prev)/norm(X) < tol, break; end
%     X_prev = X;
%   end
%
% Part 3:
%   • Set iter = k  
%   • Warn if maxIter reached without convergence
%end


function [X, iter, resHistory] = backwardGaussSeidel(A, B, tol, maxIter)

[n, m] = size(B);
if size(A,1) ~= size(A,2) || size(A,1) ~= n
    error('A must be square and match the number of rows in B.');
end

X = ones(n, m);           % initial guess
X_prev = X;                % store previous iterate for comparison
resHistory = zeros(maxIter+1, 1);  % to record residual norms for each step
resHistory(1) = norm(B - A*X, 'fro'); % initial residual (with zero guess)

for k = 1:maxIter
    for i = n:-1:1 % Backward sweep: last row to first
        if abs(A(i,i)) < eps
            error('Zero diagonal element encountered at row %d.', i);
        end
        X(i,:) = (B(i,:) - A(i,1:i-1)*X(1:i-1,:) - A(i,i+1:n)*X(i+1:n,:)) / A(i,i);
    end
    resHistory(k+1) = norm(B - A*X, 'fro');
    % Convergence check
    if norm(X - X_prev, 'fro') / max(norm(X, 'fro'), eps) < tol
        break;
    end
    X_prev = X;
end

iter = k;
if iter == maxIter && norm(X - X_prev, 'fro') / max(norm(X, 'fro'), eps) >= tol
    warning('Backward Gauss-Seidel did not converge within maxIter.');
end
resHistory = resHistory(1:iter+1); % Trim unused entries
end
