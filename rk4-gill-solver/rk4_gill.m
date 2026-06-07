% Artem Poddierohin
function [t, Y] = rk4_gill(F, t0, Y0, h, N)
% Gill's 4th-order Runge-Kutta method for a system of 2 ODE
%
% Solves the IVP:
%   Y'(t) = F(t, Y),   Y(t0) = Y0,
% where Y is a 2x1 vector (two equations)
%
% Input:
%   F - function handle, F(t, Y) returns a 2x1 vector
%   t0 - initial time
%   Y0 - initial state [y1_0; y2_0] (2x1 or 1x2)
%   h - step size 
%   N - number of steps
%
% Output:
%   t - (N+1)x1 vector of time points t_n = t0 + n*h
%   Y - (N+1)x2 matrix; row n+1 is approximation at time t(n)
%
% Method (Gill RK4):
%   K1 = h * F(tn, Yn)
%   K2 = h * F(tn + h/2, Yn + (1/2)*K1)
%   K3 = h * F(tn + h/2, Yn + ( (sqrt(2)-1)/2 )*K1 + (1 - sqrt(2)/2)*K2 )
%   K4 = h * F(tn + h, Yn - (sqrt(2)/2)*K2 + (1 + sqrt(2)/2)*K3 )
%   Yn+1 = Yn + (1/6) * ( K1 + (2-sqrt(2))*K2 + (2+sqrt(2))*K3 + K4 )

t = (t0 + (0:N)' * h);   % (N+1)x1
Y = zeros(N+1, 2);       % store as rows
Y(1, :) = Y0.';          % initial row

s2 = sqrt(2);

for n = 1:N
    tn = t(n);
    Yn = Y(n, :).';  % current state as 2x1 column

    % Stage 1
    K1 = h * F(tn, Yn);

    % Stage 2
    K2 = h * F(tn + 0.5*h, Yn + 0.5*K1);

    % Stage 3 (Gill coefficients)
    K3 = h * F(tn + 0.5*h, Yn + 0.5*(-1 + s2)*K1 + (1 - 0.5*s2)*K2);

    % Stage 4 (Gill coefficients)
    K4 = h * F(tn + h, Yn - 0.5*s2*K2 + (1 + 0.5*s2)*K3);

    % Combine stages (Gill weights)
    Yn1 = Yn + (1/6) * (K1 + (2 - s2)*K2 + (2 + s2)*K3 + K4);

    % Save result
    Y(n+1, :) = Yn1.';
end
end