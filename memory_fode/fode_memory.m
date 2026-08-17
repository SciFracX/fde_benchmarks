function fode_memory(N)
%FODE_MEMORY Run one MATLAB PECE solve for the FODE RSS benchmark.
% The benchmark is the 10-state system D_t^0.8 y = A*y on [0, 20], where A
% has -1.5 on the diagonal and 0.25 on its first off-diagonals. The initial
% state is [1.0; 0.9; ...; 0.1].
%
% Usage from MATLAB:
%   fode_memory(1280)
%   fode_memory('baseline')
%
% Runtime is measured around fde_pi12_pc only.  The parent Julia script
% measures peak RSS for this complete MATLAB process using /usr/bin/time.

    is_baseline = nargin == 1 && ...
        (ischar(N) || (isstring(N) && isscalar(N))) && ...
        strcmpi(string(N), "baseline");

    if ~is_baseline && (nargin ~= 1 || ~isnumeric(N) || ~isscalar(N) || ...
            ~isfinite(N) || N <= 0 || N ~= floor(N))
        error('fode_memory:InvalidN', 'N must be a positive integer scalar.');
    end

    this_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(this_dir);
    solver_dir = fullfile(project_dir, 'benchmarks', 'benchmark_fode', ...
        'benchmark_scipt_MATLAB');
    addpath(solver_dir);
    if isempty(which('fde_pi12_pc'))
        error('fode_memory:MissingSolver', ...
            'Cannot find fde_pi12_pc.m under %s.', solver_dir);
    end

    if is_baseline
        % MATLAB itself and the solver path are initialized, but no numerical
        % solve is performed.  This is the full-process RSS reference level.
        fprintf('FODE_MEMORY_RESULT,0,0,0\n');
        pause(1.0);
        return;
    end

    system_size = 10;
    alpha = 0.8 .* ones(system_size, 1);
    t0 = 0.0;
    T = 20.0;
    y0 = linspace(1.0, 0.1, system_size).';
    system_matrix = -1.5 .* eye(system_size) + ...
        0.25 .* diag(ones(system_size - 1, 1), 1) + ...
        0.25 .* diag(ones(system_size - 1, 1), -1);
    f_fun = @(t, y) system_matrix * y; %#ok<NASGU,INUSD>

    % Warm the MATLAB solver path before the reported runtime measurement.
    h = (T - t0) / N;
    warmup_steps = min(N, 64);
    warmup_T = t0 + warmup_steps * h;
    [warmup_t, warmup_y] = fde_pi12_pc( ...
        alpha, f_fun, t0, warmup_T, y0, h); %#ok<ASGLU>
    clear warmup_t warmup_y;

    timer_id = tic;
    [t, y] = fde_pi12_pc(alpha, f_fun, t0, T, y0, h); %#ok<ASGLU>
    runtime_seconds = toc(timer_id);

    if isempty(y) || size(y, 1) ~= system_size || ...
            size(y, 2) ~= N + 1 || any(~isfinite(y(:, end)), 'all')
        error('fode_memory:InvalidSolution', ...
            'MATLAB PECE returned a solution with invalid size or values.');
    end
    if numel(t) ~= N + 1 || abs(t(end) - T) > 64 * eps(T)
        error('fode_memory:InvalidTimeGrid', ...
            'MATLAB PECE did not return the requested fixed time grid.');
    end

    % Keep t and y alive until the machine-readable result has been emitted.
    fprintf('FODE_MEMORY_RESULT,%d,%.17g,%.17g\n', N, h, runtime_seconds);

    % R2025b can finish very short batch jobs while its asynchronous license
    % logger is still initializing.  This delay is outside the reported solver
    % runtime and reduces that startup/shutdown race without affecting RSS.
    pause(1.0);
end
