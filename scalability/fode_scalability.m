function fode_scalability(dimension, final_time, number_of_steps)
%FODE_SCALABILITY Run one MATLAB PECE scalability benchmark.
%
% Usage:
%   fode_scalability(10, 20, 5120)
%   fode_scalability('baseline')
%
% The benchmark is D_t^0.8 y = A*y, where A has -1.5 on the diagonal and
% 0.25 on its first off-diagonals. Runtime covers fde_pi12_pc only; peak RSS
% is measured for the complete process by the parent Julia script.

    is_baseline = nargin == 1 && ...
        (ischar(dimension) || (isstring(dimension) && isscalar(dimension))) && ...
        strcmpi(string(dimension), "baseline");

    if is_baseline
        initialize_solver_path();
        fprintf('SCALABILITY_RESULT,0,0,0,0,0\n');
        pause(1.0);
        return;
    end

    if nargin ~= 3
        error('fode_scalability:InvalidInput', ...
            'Usage: fode_scalability(D, T, N) or fode_scalability(''baseline'').');
    end
    validate_positive_integer(dimension, 'D');
    validate_positive_integer(number_of_steps, 'N');
    if ~isnumeric(final_time) || ~isscalar(final_time) || ...
            ~isfinite(final_time) || final_time <= 0
        error('fode_scalability:InvalidT', 'T must be a positive finite scalar.');
    end

    initialize_solver_path();

    alpha = 0.8 .* ones(dimension, 1);
    t0 = 0.0;
    y0 = linspace(1.0, 0.1, dimension).';
    system_matrix = -1.5 .* eye(dimension) + ...
        0.25 .* diag(ones(dimension - 1, 1), 1) + ...
        0.25 .* diag(ones(dimension - 1, 1), -1);
    f_fun = @(t, y) system_matrix * y; %#ok<NASGU,INUSD>
    h = (final_time - t0) / number_of_steps;

    warmup_steps = min(number_of_steps, 64);
    warmup_final_time = t0 + warmup_steps * h;
    [warmup_t, warmup_y] = fde_pi12_pc( ...
        alpha, f_fun, t0, warmup_final_time, y0, h); %#ok<ASGLU>
    clear warmup_t warmup_y;

    timer_id = tic;
    [t, y] = fde_pi12_pc( ...
        alpha, f_fun, t0, final_time, y0, h); %#ok<ASGLU>
    runtime_seconds = toc(timer_id);

    if isempty(y) || size(y, 1) ~= dimension || ...
            size(y, 2) ~= number_of_steps + 1 || ...
            any(~isfinite(y(:, end)), 'all')
        error('fode_scalability:InvalidSolution', ...
            'MATLAB PECE returned a solution with invalid size or values.');
    end
    if numel(t) ~= number_of_steps + 1 || ...
            abs(t(end) - final_time) > 64 * eps(final_time)
        error('fode_scalability:InvalidTimeGrid', ...
            'MATLAB PECE did not return the requested fixed time grid.');
    end

    fprintf('SCALABILITY_RESULT,%d,%.17g,%d,%.17g,%.17g\n', ...
        dimension, final_time, number_of_steps, h, runtime_seconds);

    % Avoid the R2025b asynchronous LicenseLogger startup/shutdown race.
    pause(1.0);
end

function initialize_solver_path()
    this_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(this_dir);
    solver_dir = fullfile(project_dir, 'benchmarks', 'benchmark_fode', ...
        'benchmark_scipt_MATLAB');
    addpath(solver_dir);
    if isempty(which('fde_pi12_pc'))
        error('fode_scalability:MissingSolver', ...
            'Cannot find fde_pi12_pc.m under %s.', solver_dir);
    end
end

function validate_positive_integer(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
            value <= 0 || value ~= floor(value)
        error('fode_scalability:InvalidInteger', ...
            '%s must be a positive integer scalar.', name);
    end
end
