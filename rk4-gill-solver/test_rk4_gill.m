% Artem Poddierohin
function test_rk4_gill()
% Testing Gill's RK4 method for systems of 2 ODE
%
% A) 5 systems with analytic (exact) solutions
% B) 5 systems with numerical reference (ode45)
%
% For each example:
%   - table of computed, reference values and absolute and relative errors
%   - one figure with y1(t), y2(t) using the SMALLEST step size

format long;

fprintf('=============================================================\n');
fprintf('Testing rk4_gill on 10 examples (A: exact, B: ode45 reference)\n\n');

h_list = [0.8 0.5 0.3 0.2 0.15 0.1 0.05 0.025];

% Example definitions

examples = {};

% A1
examples{end+1} = struct('name','A1: Exponential growth & decay', ...
    'eq','y_1'' = y_1,   y_2'' = -2 y_2', 't0',0,'T',2,'Y0',[1;2], ...
    'F',@(t,Y)[Y(1); -2*Y(2)], 'ref','exact', ...
    'exact',@(t)[exp(t), 2*exp(-2*t)] );

% A2
examples{end+1} = struct('name','A2: Constant acceleration', ...
    'eq','y_1'' = y_2,   y_2'' = 1', 't0',0,'T',3,'Y0',[0;0], ...
    'F',@(t,Y)[Y(2); 1], 'ref','exact', ...
    'exact',@(t)[0.5*t.^2, t] );

% A3
examples{end+1} = struct('name','A3: Time-dependent forcing', ...
    'eq','y_1'' = t,   y_2'' = -t', 't0',0,'T',2,'Y0',[0;1], ...
    'F',@(t,Y)[t; -t], 'ref','exact', ...
    'exact',@(t)[0.5*t.^2, 1-0.5*t.^2] );

% A4
examples{end+1} = struct('name','A4: Simple harmonic motion', ...
    'eq','y_1'' = y_2,   y_2'' = -4 y_1', 't0',0,'T',5,'Y0',[1;0], ...
    'F',@(t,Y)[Y(2); -4*Y(1)], 'ref','exact', ...
    'exact',@(t)[cos(2*t), -2*sin(2*t)] );

% A5
examples{end+1} = struct('name','A5: Nonlinear system', ...
    'eq','y_1'' = y_1^2,   y_2'' = -y_2', 't0',0,'T',0.8,'Y0',[1;2], ...
    'F',@(t,Y)[Y(1)^2; -Y(2)], 'ref','exact', ...
    'exact',@(t)[1./(1-t), 2*exp(-t)] );

% B1
examples{end+1} = struct('name','B1: Quadratic nonlinearity', ...
    'eq','y_1'' = y_2,   y_2'' = -(y_1)^2', ...
    't0',0,'T',4,'Y0',[1;0], ...
    'F',@(t,Y)[Y(2); -(Y(1)^2)], 'ref','ode45' );

% B2
examples{end+1} = struct('name','B2: Cubic damping', ...
    'eq','y_1'' = y_2,   y_2'' = -y_1 - y_2^3', ...
    't0',0,'T',10,'Y0',[1;0], ...
    'F',@(t,Y)[Y(2); -Y(1)-Y(2)^3], 'ref','ode45' );

% B3
examples{end+1} = struct('name','B3: Time forcing', ...
    'eq','y_1'' = y_2,   y_2'' = -2y_1 + sin(t)', ...
    't0',0,'T',15,'Y0',[0;0], ...
    'F',@(t,Y)[Y(2); -2*Y(1)+sin(t)], 'ref','ode45' );

% B4
examples{end+1} = struct('name','B4: Mixed nonlinear interaction', ...
    'eq','y_1'' = y_1 y_2,   y_2'' = -y_1 + y_2', ...
    't0',0,'T',2,'Y0',[1;1], ...
    'F',@(t,Y)[Y(1)*Y(2); -Y(1)+Y(2)], 'ref','ode45' );

% B5
examples{end+1} = struct('name','B5: Exponential nonlinearity', ...
    'eq','y_1'' = y_2,   y_2'' = -exp(y_1)', ...
    't0',0,'T',5,'Y0',[0;0], ...
    'F',@(t,Y)[Y(2); -exp(Y(1))], 'ref','ode45' );

% Main loop

for e = 1:length(examples)
    ex = examples{e};

    fprintf('-----------------------------------------------------\n');
    fprintf('Example %d: %s\n', e, ex.name);
    fprintf('%s\n', ex.eq);
    fprintf('t0 = %.3f, T = %.3f, Y0 = [%.6f; %.6f], reference = %s\n\n', ...
        ex.t0, ex.T, ex.Y0(1), ex.Y0(2), ex.ref);

    fprintf('%8s %8s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s\n', ...
        'h','N','y1_comp','y1_ref','y2_comp','y2_ref','absE_y1','relE_y1','absE_y2','relE_y2','absE_max','relE_max');

    for k = 1:length(h_list)
        h = h_list(k);
        N = ceil((ex.T-ex.t0)/h);
        h = (ex.T-ex.t0)/N;

        [t,Y] = rk4_gill(ex.F, ex.t0, ex.Y0, h, N);

        if strcmp(ex.ref,'exact')
            Yref = ex.exact(t(:));
        else
            opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
            sol = ode45(@(tt,YY)ex.F(tt,YY),[ex.t0 ex.T],ex.Y0,opts);
            Yref = deval(sol,t).';
        end

        err = Y-Yref;

        absE1 = max(abs(err(:,1)));
        absE2 = max(abs(err(:,2)));
        absEmax = max(absE1,absE2);

        den1 = max(abs(Yref(:,1)));
        den2 = max(abs(Yref(:,2)));
        denMax = max(max(abs(Yref)));

        if den1 > 0
            relE1 = absE1 / den1;
        else
            relE1 = NaN;
        end

        if den2 > 0
            relE2 = absE2 / den2;
        else
            relE2 = NaN;
        end

        if denMax > 0
            relEmax = absEmax / denMax;
        else
            relEmax = NaN;
        end

        y1c = Y(end,1);
        y1r = Yref(end,1);
        y2c = Y(end,2);
        y2r = Yref(end,2);

        fprintf('%8.4f %8d %14.6e %14.6e %14.6e %14.6e %14.3e %14.3e %14.3e %14.3e %14.3e %14.3e\n', ...
            h,N,y1c,y1r,y2c,y2r,absE1,relE1,absE2,relE2,absEmax,relEmax);
    end

    % plot with smallest h
    h = h_list(end);
    N = ceil((ex.T-ex.t0)/h);
    h = (ex.T-ex.t0)/N;

    [t,Y] = rk4_gill(ex.F, ex.t0, ex.Y0, h, N);

    if strcmp(ex.ref,'exact')
        Yref = ex.exact(t(:));
    else
        opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
        sol = ode45(@(tt,YY)ex.F(tt,YY),[ex.t0 ex.T],ex.Y0,opts);
        Yref = deval(sol,t).';
    end

    figure('Name',ex.name);

    subplot(2,1,1)
    plot(t,Y(:,1),'b-',t,Yref(:,1),'r--')
    grid on
    ylabel('y_1(t)')
    title(ex.eq)
    legend('Gill RK4','Reference')

    subplot(2,1,2)
    plot(t,Y(:,2),'b-',t,Yref(:,2),'r--')
    grid on
    xlabel('t')
    ylabel('y_2(t)')
    legend('Gill RK4','Reference')
end

fprintf('\nAll examples completed.\n');
end
