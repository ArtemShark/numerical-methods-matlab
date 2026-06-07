% Artem Poddierohin
%function [X, iter, resHistory] = gaussSeidel(A, B, tol, maxIter)
% gaussSeidel   Solve A*X = B using the forward Gauss–Seidel iteration
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
%   • Check that A is square and size(B,1)==n  
%   • Pre-allocate X = zeros(n,m) and resHistory(1) = norm(B–A*X)  
%
% Part 2: Iterative loop
%   for k = 1:maxIter Each iteration k is one full pass over all rows
%     for i = 1:n                      % forward sweep
%       X(i,:) = (B(i,:) 
%                – A(i,1:i-1)*X(1:i-1,:) multiplies the newly updated entries above the diagonal
%                – A(i,i+1:n)*X(i+1:n,:)) / A(i,i); uses the old entries from the previous sweep
%     end
%     resHistory(k+1) = norm(B–A*X);  (track how fast ||B - AX^(k)|| is shrinking)
%     if norm(X–X_prev)/norm(X) < tol, break; end
%     (Compute the relative change between the new and old iterates, 
%     If this falls below tol, we break out of the loop early.)
%     X_prev = X; (stores the current iterate so the next pass can compare it again.)
%   end
%
% Part 3: 
%   • Set iter = k  (Return actual iteration count)
%   • If iter == maxIter and not converged, issue a warning
%   (reached maxIter without satisfying the tolerance)
%end


function [X, iter, resHistory] = gaussSeidel(A, B, tol, maxIter)

[n, m] = size(B);
if size(A,1) ~= size(A,2) || size(A,1) ~= n
    error('A must be square and match the number of rows in B.');
end

X = zeros(n, m);             % initial guess
X_prev = X;                  % store previous iterate for comparison
resHistory = zeros(maxIter+1, 1);  % to record residual norms for each step
resHistory(1) = norm(B - A*X, 'fro');  % initial residual (with zero guess)

for k = 1:maxIter
    for i = 1:n
        if abs(A(i,i)) < eps
            error('Zero diagonal element encountered at row %d.', i);
        end
        % Update row i using the most recent values for X
        X(i,:) = (B(i,:) - A(i,1:i-1)*X(1:i-1,:) - A(i,i+1:n)*X(i+1:n,:)) / A(i,i);
    end
    resHistory(k+1) = norm(B - A*X, 'fro');
    % Check for convergence (relative change)
    if norm(X - X_prev, 'fro') / max(norm(X, 'fro'), eps) < tol
        break;
    end
    X_prev = X;
end

iter = k;
if iter == maxIter && norm(X - X_prev, 'fro') / max(norm(X, 'fro'), eps) >= tol
    warning('Gauss-Seidel did not converge within maxIter.');
end
resHistory = resHistory(1:iter+1); % Trim unused entries
end