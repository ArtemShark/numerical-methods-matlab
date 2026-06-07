% Artem Poddierohin
function Integral = approx_trapez_rule(f, Nr, Ntheta)
% Approximate integral over the unit disk x^2 + y^2 <= 1
%
% Integral = approx_trapez_rule(f, Nr, Ntheta)
%
% Input:
%    f - handle to function of two variables
%    Nr - number of subintervals in r (0 <= r <= 1)
%    Ntheta - number of subintervals in θ (0 <= θ <= 2π)
% Output:
%    Integral - approximate value of ∫∫{x^2+y^2<=1} f(x,y) dx dy
%
% Method:
%    We compute ∫∫{x^2+y^2<=1} f(x,y) dx dy using polar coordinates:
%           x = r*cos(θ),  y = r*sin(θ),
%           0 <= r <= 1,      0 <= θ <= 2π,
%           dx dy = r dr dθ
%    Then we apply the 2D composite trapezoidal rule in r and θ

% Integration intervals in polar coordinates
r_min = 0;
r_max = 1;
t_min = 0;
t_max = 2*pi;

% Step sizes
hr = (r_max - r_min) / Nr;
ht = (t_max - t_min) / Ntheta;

% Grid points
r = r_min + (0:Nr) * hr;        % r_0, r_1, ..., r_Nr
t = t_min + (0:Ntheta) * ht;    % θ_0, ..., θ_Ntheta

% Trapezoidal weights in r: 1/2 on the ends, 1 inside
wr = ones(1, Nr+1);
wr(1)   = 0.5;
wr(end) = 0.5;

% Trapezoidal weights in θ: 1/2 on the ends, 1 inside
wt = ones(1, Ntheta+1);
wt(1)   = 0.5;
wt(end) = 0.5;

% Double sum for the 2D composite trapezoidal rule
Integral_sum = 0;

for i = 1:(Nr+1)
    ri = r(i);  % current radius

    for j = 1:(Ntheta+1)
        tj = t(j);  % current angle

        % Transform to Cartesian coordinates
        x = ri * cos(tj);
        y = ri * sin(tj);

        % Value of f(x,y)
        fx = f(x, y);

        % Integrand in polar coordinates: f(x,y) * r
        G = fx * ri;

        % Weighted contribution with trapezoidal weights
        Integral_sum = Integral_sum + wr(i) * wt(j) * G;
    end
end

% Multiply by step sizes hr and ht
Integral = hr * ht * Integral_sum;
end