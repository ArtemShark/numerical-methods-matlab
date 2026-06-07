% Artem Poddierohin
function test_GS_methods()
tol = 1e-10; maxIter = 1000;
rng(101); % Reproducible random numbers

% ========== SMALL WELL-CONDITIONED EXAMPLES ==========
dims = [6 6 7 7 8];
fprintf('======= SMALL WELL-CONDITIONED EXAMPLES =======\n');
for k = 1:length(dims)
    n = dims(k);
    A = randn(n); 
    A = A + diag(sum(abs(A),2));
    X_true = randi([-5 5], n, 1);  % Different X_true for each test
    B = A * X_true;

    [X_gs, iter_gs, ~] = gaussSeidel(A, B, tol, maxIter);
    [X_bgs, iter_bgs, ~] = backwardGaussSeidel(A, B, tol, maxIter);

    fprintf('\n------ Small Example %d: %dx%d ------\n', k, n, n);
    printMatrix('Matrix A:', A);
    printVector('Right-hand side B:', B);
    printVector('Expected solution X_true:', X_true);
    printVector('Gauss-Seidel computed X:', X_gs);
    printVector('Backward GS computed X:', X_bgs);
    fprintf('GS iters: %d, BGS iters: %d\n', iter_gs, iter_bgs);
    fprintf('GS rel.err: %.2e, BGS rel.err: %.2e\n', ...
        norm(X_gs-X_true)/norm(X_true), norm(X_bgs-X_true)/norm(X_true));
end

% ========== SMALL FAILURE CASES ==========
fprintf('\n======= SMALL FAILURE CASES =======\n');
for k = 1:2
    n = 6;
    A = randn(n); % Not diagonally dominant
    if k == 2
        A(2,2) = 0; % Zero diagonal for 2nd fail case
    end
    X_true = randi([-3 3], n, 1);
    B = A * X_true;

    fprintf('\nFailure Case %d: Random 6x6, non-dominant or singular\n', k);
    printMatrix('Matrix A:', A);
    printVector('Right-hand side B:', B);
    printVector('Expected solution X_true:', X_true);
    try
        [X_gs, iter_gs, ~] = gaussSeidel(A, B, tol, maxIter);
        printVector('Gauss-Seidel computed X:', X_gs);
        fprintf('GS iters: %d\n', iter_gs);
        fprintf('GS rel.err: %.2e\n', norm(X_gs-X_true)/norm(X_true));
    catch ME
        fprintf('Gauss-Seidel failed: %s\n', ME.message);
    end
    try
        [X_bgs, iter_bgs, ~] = backwardGaussSeidel(A, B, tol, maxIter);
        printVector('Backward GS computed X:', X_bgs);
        fprintf('BGS iters: %d\n', iter_bgs);
        fprintf('BGS rel.err: %.2e\n', norm(X_bgs-X_true)/norm(X_true));
    catch ME
        fprintf('Backward GS failed: %s\n', ME.message);
    end
end

% ========== LARGE RANDOM CASES (SUMMARY TABLE) ==========
fprintf('\n======= LARGE RANDOM CASES (SUMMARY TABLE) =======\n');
fprintf([' n  cond(A)    rho_GS  rho_BGS  iter_GS  iter_BGS   time_GS(s)   time_BGS(s)', ...
    '   relErr_GS   relErr_BGS  fwdErr_GS  fwdErr_BGS  backErr_GS  backErr_BGS\n']);
for n = 100:100:800
    A = randn(n) + n*eye(n);  % Well-conditioned
    X_true = randn(n,1);
    B = A * X_true;
    
    t1 = tic;
    [X_gs, iter_gs, ~] = gaussSeidel(A, B, tol, maxIter);
    t_gs = toc(t1);

    t1 = tic;
    [X_bgs, iter_bgs, ~] = backwardGaussSeidel(A, B, tol, maxIter);
    t_bgs = toc(t1);

    condA = cond(A);
    % Iteration matrices
    D = diag(diag(A)); L = tril(A,-1); U = triu(A,1);
    T_gs = -(D+L)\U; 
    T_bgs = -(D+U)\L;
    rho_gs = max(abs(eig(T_gs)));
    rho_bgs = max(abs(eig(T_bgs)));

    err_gs = norm(X_gs-X_true)/norm(X_true);
    err_bgs = norm(X_bgs-X_true)/norm(X_true);
    fwd_gs = err_gs / condA;
    fwd_bgs = err_bgs / condA;
    back_gs = norm(B - A*X_gs, 'fro')/(norm(A, 'fro')*norm(X_gs, 'fro'));
    back_bgs = norm(B - A*X_bgs, 'fro')/(norm(A, 'fro')*norm(X_bgs, 'fro'));

    fprintf('%3d %9.2e %8.3f %8.3f %8d %9d %11.4f %11.4f %11.2e %11.2e %11.2e %11.2e %11.2e %11.2e\n', ...
        n, condA, rho_gs, rho_bgs, iter_gs, iter_bgs, t_gs, t_bgs, ...
        err_gs, err_bgs, fwd_gs, fwd_bgs, back_gs, back_bgs);
end


fprintf('\nAll test cases complete.\n');
end

% ====== Helper functions =====

function printMatrix(titleStr, A)
fprintf('%s\n', titleStr);
for r = 1:size(A,1)
    fprintf('%10.4f', A(r,:));
    fprintf('\n');
end
end

function printVector(titleStr, v)
fprintf('%s\n', titleStr);
fprintf('%10.4f', v);
fprintf('\n');
end