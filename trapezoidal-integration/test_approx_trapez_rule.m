% Artem Poddierohin 
function test_approx_trapez_rule()
% Testing for approx_trapez_rule.m
%
% For each test function f(x,y):
%   1) we compute the exact value of the integral over the unit disk
%   2) we approximate the integral using approx_trapez_rule,
%   3) we print absolute and relative errors in a table
%
% Region: D = { (x,y): x^2 + y^2 <= 1 }

format long;

% Mesh sizes for testing
Nr_vals     = [ 4  8  16  32   8  16  12  24  20  30];
Ntheta_vals = [ 4  8  16  32  16   8  24  12  30  20];

% Common polar limits
r_min = 0;
r_max = 1;
t_min = 0;
t_max = 2*pi;

% Table formatting
headerFmt = '%6s  %6s  %20s  %12s  %12s\n';
rowFmt    = '%6d  %6d  %20.15f  %12.3e  %12.3e\n';

% Example 1: f(x,y) = 1

% Inner integral in r: ∫0^1 r dr = (r^2)/2 |0^1
I_r1 = (r_max^2 - r_min^2) / 2;
% Outer integral in theta: ∫0^2π dθ = 2π
I_t1 = (t_max - t_min);
% Exact integral: Integral = (2π) * (1/2) = π
Integral_exact1 = I_t1 * I_r1;

fprintf('=== Example 1: f(x,y) = 1 ===\n');
fprintf('Exact value: Integral_exact1 = %.15f\n', Integral_exact1);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f1 = @(x,y) 1; 
    Integral_num = approx_trapez_rule(f1, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact1);
    relErr = absErr / abs(Integral_exact1);
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 2: f(x,y) = x^2 + y^2

% In polar: x^2 + y^2 = r^2.
% Integral = ∫0^2π ∫0^1 r^2 * r dr dθ = ∫0^2π ∫0^1 r^3 dr dθ
% Inner: ∫0^1 r^3 dr = (r^4)/4 |0^1
I_r2 = (r_max^4 - r_min^4) / 4;
% Outer: ∫0^2π dθ = 2π
I_t2 = (t_max - t_min);
% Exact integral: Integral = 2π * (1/4) = π/2
Integral_exact2 = I_t2 * I_r2;

fprintf('=== Example 2: f(x,y) = x^2 + y^2 ===\n');
fprintf('Exact value: Integral_exact2 = %.15f\n', Integral_exact2);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f2 = @(x,y) x.^2 + y.^2;
    Integral_num = approx_trapez_rule(f2, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact2);
    relErr = absErr / abs(Integral_exact2);
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 3: f(x,y) = sin(x^2 + y^2)

% In polar: x^2 + y^2 = r^2  => f = sin(r^2)
% Integral = ∫0^2π ∫0^1 sin(r^2) * r dr dθ
% u = r^2, du = 2r dr, so r dr = du/2
% Inner:
%   ∫0^1 sin(r^2) r dr = (1/2) ∫0^1 sin(u) du = (1/2)[ -cos(u) ]0^1 = (1/2)(1 - cos(1))
% Outer: ∫0^2π dθ = 2π => Integral = 2π * (1/2)(1 - cos(1)) = π(1 - cos(1))
u_min3 = r_min^2;     % 0
u_max3 = r_max^2;     % 1
I_r3_inner = 0.5 * (-cos(u_max3) + cos(u_min3));  % (1/2)*∫ sin(u) du
I_t3 = (t_max - t_min);                           % 2π
Integral_exact3 = I_t3 * I_r3_inner;

fprintf('=== Example 3: f(x,y) = sin(x^2 + y^2) ===\n');
fprintf('Exact value: Integral_exact3 = %.15f\n', Integral_exact3);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f3 = @(x,y) sin(x.^2 + y.^2);
    Integral_num = approx_trapez_rule(f3, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact3);
    relErr = absErr / abs(Integral_exact3);
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 4: f(x,y) = exp(-(x^2 + y^2))

% In polar: f = exp(-r^2).
% Integral = ∫0^2π ∫0^1 e^{-r^2} * r dr dθ
% u = r^2, du = 2r dr, r dr = du/2.
% Inner:
%   ∫0^1 e^{-r^2} r dr = (1/2) ∫0^1 e^{-u} du = (1/2)[ -e^{-u} ]0^1 = (1/2)(1 - e^{-1})
% Outer: ∫0^2π dθ = 2π => Integral = 2π * (1/2)(1 - e^{-1}) = π(1 - e^{-1})
u_min4 = r_min^2;     % 0
u_max4 = r_max^2;     % 1
I_r4_inner = 0.5 * (-exp(-u_max4) + exp(-u_min4));  % (1/2)*∫ e^{-u} du
I_t4 = (t_max - t_min);                             % 2π
Integral_exact4 = I_t4 * I_r4_inner;

fprintf('=== Example 4: f(x,y) = exp(-(x^2 + y^2)) ===\n');
fprintf('Exact value: Integral_exact4 = %.15f\n', Integral_exact4);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f4 = @(x,y) exp(-(x.^2 + y.^2));
    Integral_num = approx_trapez_rule(f4, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact4);
    relErr = absErr / abs(Integral_exact4);
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 5: f(x,y) = x^2 - y^2 (integral = 0)

% By symmetry of the disk, positive and negative parts cancel out
Integral_exact5 = 0;

fprintf('=== Example 5: f(x,y) = x^2 - y^2 ===\n');
fprintf('Exact value: Integral_exact5 = %.15f\n', Integral_exact5);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f5 = @(x,y) x.^2 - y.^2;
    Integral_num = approx_trapez_rule(f5, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact5);
    % Relative error not meaningful for exact = 0
    relErr = NaN;
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 6: f(x,y) = (x^2 + y^2)^5  

k = 5; % (x^2 + y^2)^k -> degree 2k
Integral_exact6 = pi / (k + 1); % = pi/6

fprintf('=== Example 6: f(x,y) = (x^2 + y^2)^5 ===\n');
fprintf('Exact value: Integral_exact6 = %.15f\n', Integral_exact6);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for idx = 1:length(Nr_vals)
    Nr = Nr_vals(idx);
    Nt = Ntheta_vals(idx);
    f6 = @(x,y) (x.^2 + y.^2).^5;
    Integral_num = approx_trapez_rule(f6, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact6);
    relErr = absErr / abs(Integral_exact6);
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 7: f(x,y) = cos(x^2 + y^2)

% In polar: f = cos(r^2)
% Integral = ∫0^2π ∫0^1 cos(r^2) * r dr dθ
% u = r^2, du = 2r dr => r dr = du/2
% Inner: (1/2)∫_0^1 cos(u) du = (1/2)[sin(u)]_0^1 = (1/2)sin(1)
% Outer: ∫0^2π dθ = 2π => Integral = 2π * (1/2)sin(1) = π sin(1)
u_min7 = r_min^2;
u_max7 = r_max^2;
I_r7_inner = 0.5 * (sin(u_max7) - sin(u_min7));
I_t7 = (t_max - t_min);
Integral_exact7 = I_t7 * I_r7_inner;

fprintf('=== Example 7: f(x,y) = cos(x^2 + y^2) ===\n');
fprintf('Exact value: Integral_exact7 = %.15f\n', Integral_exact7);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f7 = @(x,y) cos(x.^2 + y.^2);
    Integral_num = approx_trapez_rule(f7, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact7);
    relErr = absErr / abs(Integral_exact7);
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 8: f(x,y) = x^2 * y^2

% In polar: x^2 y^2 = r^4 cos^2(θ) sin^2(θ)
% Integrand in polar: r^5 cos^2(θ) sin^2(θ)
% Radial: ∫0^1 r^5 dr = 1/6
% Angular: cos^2(θ) sin^2(θ) = (1/4) sin^2(2θ)
%   ∫0^2π sin^2(2θ) dθ = (average 1/2) * 2π = π
%   ∫0^2π cos^2(θ) sin^2(θ) dθ = (1/4)*π = π/4
% Total: (1/6)*(π/4) = π/24
I_r8 = (r_max^6 - r_min^6) / 6;
I_t8 = pi/4;
Integral_exact8 = I_t8 * I_r8;

fprintf('=== Example 8: f(x,y) = x^2 * y^2 ===\n');
fprintf('Exact value: Integral_exact8 = %.15f\n', Integral_exact8);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f8 = @(x,y) x.^2 .* y.^2;
    Integral_num = approx_trapez_rule(f8, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact8);
    relErr = absErr / abs(Integral_exact8);
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 9: f(x,y) = x^5  (integral = 0)

% The disk is symmetric in x -> -x, while x^5 is odd in x,
% so the integral over the disk is 0 by symmetry
Integral_exact9 = 0;

fprintf('=== Example 9: f(x,y) = x^5 ===\n');
fprintf('Exact value: Integral_exact9 = %.15f\n', Integral_exact9);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f9 = @(x,y) x.^5;
    Integral_num = approx_trapez_rule(f9, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact9);
    relErr = NaN;  % relative error not meaningful when exact = 0
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

% Example 10: f(x,y) = ln(1 + x^2 + y^2)

% In polar: f = ln(1 + r^2)
% Integral = ∫0^2π ∫0^1 ln(1 + r^2) * r dr dθ
% u = 1 + r^2, du = 2 r dr -> r dr = du/2
% Inner: (1/2) ∫1^2 ln(u) du = (1/2)[u ln u - u]1^2 = (1/2)(2 ln 2 - 2 + 1) = (1/2)(2 ln 2 - 1)
% Outer: ∫0^2π dθ = 2π -> Integral = π (2 ln 2 - 1)

u_min10 = 1 + r_min^2;   % 1
u_max10 = 1 + r_max^2;   % 2
I_r10 = 0.5 * ((u_max10*log(u_max10) - u_max10) - (u_min10*log(u_min10) - u_min10));
I_t10 = (t_max - t_min);
Integral_exact10 = I_t10 * I_r10;   % π(2 ln 2 - 1)

fprintf('=== Example 10: f(x,y) = ln(1 + x^2 + y^2) ===\n');
fprintf('Exact value: Integral_exact10 = %.15f\n', Integral_exact10);
fprintf(headerFmt, 'Nr', 'Ntheta', 'approx value', 'abs.err', 'rel.err');

for k = 1:length(Nr_vals)
    Nr = Nr_vals(k);
    Nt = Ntheta_vals(k);
    f10 = @(x,y) log(1 + x.^2 + y.^2); 
    Integral_num = approx_trapez_rule(f10, Nr, Nt);
    absErr = abs(Integral_num - Integral_exact10);
    relErr = absErr / abs(Integral_exact10);
    fprintf(rowFmt, Nr, Nt, Integral_num, absErr, relErr);
end
fprintf('\n');

fprintf('All test cases completed.\n');
end