clc
% Static steering model: input magnetic induction intensity B0 (Tesla) and tilted angle of magnetic field alpha (deg),
% output static deflection angle of magnetically jammed tip theta (deg).
% For parameter sweeps, do not make alpha and B0 both vectors.

%==============================================================
% 1) Inputs
B0 = (1:1:70)'/1000; % Magnetic induction intensity in Tesla (used for physics)
% B0 = 70/1000; % Alternatively, use a scalar (Tesla)
% alpha = (1:1:89)'; % Alternatively, sweep alpha (deg)
alpha = 30; % Tilted angle of magnetic field (deg)

% For plotting/export: convert B0 (T) -> B0_mT (mT)
B0_mT = B0 .* 1000;

%==============================================================
% 2) Catheter parameters, compute Q
D1 = 4e-3; % Outer tube OD
D2 = 3.8e-3; % Outer tube ID
D3 = 1.5e-3; % Inner tube OD
D4 = 1e-3; % Inner tube ID
L = 50e-3; % Length
A = pi/4 .* (D2.^2 - D3.^2); % Cavity base area
I1 = pi/64 .* (D1.^4 - D2.^4); % Second moment of area (outer tube)
I2 = pi/64 .* (D3.^4 - D4.^4); % Second moment of area (inner tube)
G = 1e6; % Shear modulus of matrix
niu = 0.5; % Poisson's ratio of matrix
V = A .* L; % Volume of catheter cavity
Q = 3 .* pi .* (1 + niu) .* G .* (D1.^4 - D2.^4 + D3.^4 - D4.^4) ./ (32 .* L); % Custom elastic coefficient Q

%==============================================================
% 3) Magnetic particles, compute P
H = B0 .* 1000 ./ (0.4 .* pi); % Magnetic field strength in kA/m
r = 125e-6; % Particle radius
fai = 0.523; % Particle volume fraction
kai0 = 7.1; % Initial susceptibility
miu0 = 4 .* pi .* 1e-7; % Vacuum permeability
Ms = 2.156 ./ miu0; % Saturation magnetization of iron powder
M = kai0 .* Ms .* H .* 1000 ./ (kai0 .* H .* 1000 + Ms); % Frohlich–Kennelly law
m = M .* (4/3) .* pi .* r.^3; % Magnetic dipole moment of one particle
miu = 0.05; % Friction coefficient between particles
P = miu0 .* fai .* pi ./ 96 .* (D2.^2 - D3.^2) .* (kai0 .* Ms .* B0 ./ (miu0 .* Ms + kai0 .* B0)).^2; % Custom magnetic coefficient P
Q_P_ratio = Q ./ P; % Q/P ratio

%==============================================================
% 4) Numerical solution for theta
% Solve elementwise: P(i)*sin(2(alpha(i)-theta(i))) = Q(i)*sin(theta(i))
% alpha, P, Q can be scalars or vectors; alpha/theta are in degrees; solve in radians internally.

% 4.1 Normalize to column vectors
makecol = @(x) reshape(x, [], 1);
alpha = makecol(alpha);
P_in = P; Q_in = Q; % backup
if ~isscalar(P), P = makecol(P); end
if ~isscalar(Q), Q = makecol(Q); end

% 4.2 Target length
lens = [numel(alpha), numel(P), numel(Q)];
n = max(lens);

% 4.3 Broadcast scalars
if isscalar(alpha) && n > 1, alpha = repmat(alpha, n, 1); end
if isscalar(P) && n > 1, P = repmat(P, n, 1); end
if isscalar(Q) && n > 1, Q = repmat(Q, n, 1); end

% 4.4 Size consistency check
lens = [numel(alpha), numel(P), numel(Q)];
if any(lens ~= n)
fprintf('Size(alpha)=%s, Size(P)=%s, Size(Q)=%s\n', mat2str(size(alpha)), mat2str(size(P_in)), mat2str(size(Q_in)));
error('alpha, P, Q must be scalars or vectors of the same length (after row->column normalization).');
end

% 4.5 Solve per element (prefer fsolve, fallback to fzero)
alph = deg2rad(alpha); % radians
theta = nan(n,1); % numerical solution (deg)

% Tolerances (radians) and per-point tolerance tracking
fsolve_func_tol = 1e-12; % FunctionTolerance
fsolve_step_tol = 1e-12; % StepTolerance
default_fzero_tol_rad = 1e-10; % Conservative angle tolerance if fzero is used
rad2deg_const = 180/pi;
theta_tol_rad = nan(n,1); % per-point tolerance (rad)

% First-order initial guess
den = 2 .* P + Q .* cos(alph);
mask = abs(den) < eps;
den(mask) = den(mask) + eps .* sign(den(mask) + (den(mask) == 0));
delta1 = (Q .* sin(alph)) ./ den;
theta0 = alph - delta1;
theta0 = min(max(theta0, 1e-9), alph - 1e-9); % ensure within (0, alpha)

opts = optimoptions('fsolve', ...
'Display', 'off', ...
'FunctionTolerance', fsolve_func_tol, ...
'OptimalityTolerance', 1e-12, ...
'StepTolerance', fsolve_step_tol, ...
'MaxIterations', 400);

for i = 1:n
a_i = alph(i); P_i = P(i); Q_i = Q(i);
F_i = @(th) P_i .* sin(2 .* (a_i - th)) - Q_i .* sin(th);


[th_i, ~, exitflag] = fsolve(F_i, theta0(i), opts);

used_fzero = false;
if ~(exitflag > 0 && th_i > 0 && th_i < a_i)
    % Bracket by scanning + fzero fallback
    Nscan = 200;
    lo = a_i .* 1e-9; hi = a_i .* (1 - 1e-9);
    tt = linspace(lo, hi, Nscan);
    Ft = arrayfun(F_i, tt);
    k = find(Ft(1:end-1) .* Ft(2:end) <= 0, 1, 'first');
    if ~isempty(k)
        th_i = fzero(F_i, [tt(k), tt(k+1)]);
    else
        th_i = fzero(F_i, theta0(i));
    end
    used_fzero = true;
end

if ~(th_i > 0 && th_i < a_i)
    error('Failed to solve at index %d (alpha=%.6g deg, P=%.6g, Q=%.6g).', i, rad2deg(a_i), P_i, Q_i);
end

% Record per-point tolerance
if used_fzero
    theta_tol_rad(i) = default_fzero_tol_rad;
else
    theta_tol_rad(i) = max(fsolve_step_tol, fsolve_func_tol);
end

theta(i) = rad2deg(th_i);
end

% Representative (conservative) tolerance in degrees
tol_deg_out = max(theta_tol_rad) .* rad2deg_const;

%==============================================================
% 5) Plot and export variables
% Notes:
% - If B0 is a vector (and alpha is scalar), x_out uses B0_mT (mT)
% - If alpha is a vector (and B0 is scalar), x_out uses alpha (deg)

% ================== Add-on: Comparison along B0 or alpha (export-friendly variables) ==================
% Workspace prerequisites: alpha, B0, B0_mT, P, Q, theta (theta in degrees)
toCol = @(x) reshape(x, [], 1);
isConstVec = @(x) numel(x) > 1 && all(x(:) == x(1)); % Treat an all-equal vector as a constant

% Safe length-align helper (no ternary operator)
function y = toLen_safe(x, n)
if isscalar(x)
y = repmat(x, n, 1);
else
y = reshape(x, [], 1);
if numel(y) ~= n
% If mismatch but x is a constant vector, compress to scalar and expand
if numel(y) > 1 && all(y == y(1))
y = repmat(y(1), n, 1);
else
error('toLen_safe: length mismatch. Expected %d, got %d.', n, numel(y));
end
end
end
end

% Normalize shapes (do not alter original variables)
a = toCol(alpha);
bT = toCol(B0); % B0 in Tesla (for physics)
b_mT = toCol(B0_mT); % B0 in mT (for plotting/export)
thNum = toCol(theta); % numerical theta (deg)

% Demote constant vectors to scalar viewpoint
alpha_is_const = isConstVec(a);
B0_is_const = isConstVec(bT);
if alpha_is_const, a = a(1); end
if B0_is_const, bT = bT(1); b_mT = b_mT(1); end

% Determine which dimension is scanned
isA_vec = numel(a) > 1;
isB_vec = numel(bT) > 1;

if ~(xor(isA_vec, isB_vec))
if ~isA_vec && ~isB_vec
fprintf('Both alpha and B0 are constants. Skip plotting. theta = %.6g deg\n', thNum);
else
warning('alpha and B0 cannot both be non-constant vectors. Please scan only one of them.');
end
% Also clear intended outputs to avoid confusion
clear x_out theta_num_out theta_approx_out abs_err_out rel_err_out
return
end

% Build arrays aligned to the scanned variable and compute first-order approximation
if isA_vec
% Scan over alpha -> x uses alpha (deg)
x = toCol(a); xlab = '\alpha (deg)';
n = numel(x);
Pv = toLen_safe(P, n);
Qv = toLen_safe(Q, n);
thv = toLen_safe(thNum, n);


alph_rad = deg2rad(x);
den = 2 .* Pv + Qv .* cos(alph_rad);
m   = abs(den) < eps; den(m) = den(m) + eps .* sign(den(m) + (den(m) == 0));
th1 = rad2deg(alph_rad - (Qv .* sin(alph_rad)) ./ den);
else
% Scan over B0 -> x uses B0 in mT
x = toCol(b_mT); xlab = 'B_0 (mT)';
n = numel(x);
Pv = toLen_safe(P, n);
Qv = toLen_safe(Q, n);
thv = toLen_safe(thNum, n);


% Use alpha as a constant (or aligned if main code already expanded it)
if isscalar(alpha) || isConstVec(alpha)
    alph_rad = deg2rad(alpha(1)) .* ones(n,1);
else
    alph_rad = deg2rad(toLen_safe(alpha, n));
end
den = 2 .* Pv + Qv .* cos(alph_rad);
m   = abs(den) < eps; den(m) = den(m) + eps .* sign(den(m) + (den(m) == 0));
th1 = rad2deg(alph_rad - (Qv .* sin(alph_rad)) ./ den);
end

% Sort for a clean curve and consistent export
[xs, idx] = sort(x);
th_num_s = thv(idx);
th1_s = th1(idx);

% Errors
abs_err = abs(th_num_s - th1_s);
rel_err = abs_err ./ max(abs(th_num_s), 1e-12);

% ================== Named outputs for export (Excel-ready) ==================
% x_out: x-axis (alpha [deg] or B0 [mT], sorted)
% theta_num_out: numerical solution theta on the same x grid [deg]
% theta_approx_out: first-order approximation on the same x grid [deg]
% abs_err_out: absolute error |theta_num − theta_approx| [deg]
% rel_err_out: relative error (unitless); multiply by 100 for percentage
% tol_deg_out: numerical-solution tolerance (deg), conservative upper bound for the sweep
x_out = xs; % alpha (deg) or B0 (mT)
theta_num_out = th_num_s; % numerical theta (deg)
theta_approx_out = th1_s; % first-order approx theta (deg)
abs_err_out = abs_err; % absolute error (deg)
rel_err_out = rel_err; % relative error (unitless)
% tol_deg_out is already computed above (scalar, deg)

% Optional: quick preview plots
figure('Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(x_out, theta_num_out, 'b-', 'LineWidth', 1.8); hold on;
plot(x_out, theta_approx_out, 'r--', 'LineWidth', 1.8);
grid on; xlabel(xlab); ylabel('\theta (deg)');
title('\theta: Numerical vs First-order Approximation');
legend('Numerical','First-order approx','Location','best');

nexttile;
yyaxis left; plot(x_out, abs_err_out, 'k-', 'LineWidth', 1.6); ylabel('Absolute error | \Delta\theta | (deg)');
yyaxis right; plot(x_out, rel_err_out .* 100, 'm--', 'LineWidth', 1.6); ylabel('Relative error (%)');
grid on; xlabel(xlab); title('Error comparison');

fprintf('Error stats: max |Δθ| = %.4g deg, median |Δθ| = %.4g deg, max relative error = %.3g %%\n', ...
max(abs_err_out), median(abs_err_out), max(rel_err_out) * 100);
fprintf('Reported numerical-solution tolerance (deg): tol_deg_out = %.3e deg\n', tol_deg_out);

% ======================================================================
% Export to Excel
% x_out is alpha[deg] or B0[mT], depending on which was scanned
T = table(x_out, theta_num_out, theta_approx_out, abs_err_out, rel_err_out, ...
repmat(tol_deg_out, numel(x_out), 1), ...
'VariableNames', {'x','theta_num_deg','theta_approx_deg','abs_err_deg','rel_err','tol_deg'});
writetable(T, 'curves.xlsx', 'Sheet', 1, 'Range', 'A1');