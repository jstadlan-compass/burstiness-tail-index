function dist = plcutoff_dist(alpha, lambda, xmin)
% PLCUTOFF_DIST  Factory for power-law with exponential cutoff distribution.
%
%   f(x) ∝ x^(-alpha) * exp(-lambda * x),   x >= xmin > 0
%
% Usage:
%   dist = plcutoff_dist(alpha, lambda, xmin);
%   dist.pdf(x)
%   dist.cdf(x)
%   dist.ppf(u)
%   dist.rnd(n)      % n samples
%   dist.mean
%   dist.var
%   dist.cv
%
% This uses the x_min version (support [xmin, ∞)) and computes everything
% via upper incomplete gamma Γ(s,z), implemented so that negative s are OK
% (except at gamma poles s = 0, -1, -2, ...).

    % ---- Basic checks ----
    if lambda <= 0
        error('lambda must be > 0.');
    end
    if xmin <= 0
        error('xmin must be > 0.');
    end

    % Shapes for incomplete gamma
    s1 = 1 - alpha;   % for Γ(1-α, ·)
    s2 = 2 - alpha;   % for Γ(2-α, ·)
    s3 = 3 - alpha;   % for Γ(3-α, ·)

    % Guard against gamma poles: s1 = 0, -1, -2, ...
    if abs(s1 - round(s1)) < 1e-12 && round(s1) <= 0
        error(['This implementation cannot handle alpha such that 1-alpha ' ...
               'is a non-positive integer (gamma pole). Choose a different alpha.']);
    end

    zmin = lambda * xmin;

    % ---- Normalizing constant and moments (at build time) ----
    % Gk = Γ(k-α, λ xmin)
    G1 = upper_incomplete_gamma(s1, zmin);   % Γ(1-α, λ xmin)
    G2 = upper_incomplete_gamma(s2, zmin);   % Γ(2-α, λ xmin)
    G3 = upper_incomplete_gamma(s3, zmin);   % Γ(3-α, λ xmin)

    if ~all(isfinite([G1, G2, G3]))
        error('Incomplete gamma evaluation failed; check parameters.');
    end

    % Normalizing constant:
    %   C = 1 / [ λ^(α-1) * Γ(1-α, λ xmin) ]
    C = 1 / (lambda^(alpha - 1) * G1);

    % Moments:
    %   E[X^n] = λ^(-n) * Γ(n+1-α, λxmin) / Γ(1-α, λxmin)
    meanX = (lambda^-1) * (G2 / G1);
    EX2   = (lambda^-2) * (G3 / G1);
    varX  = EX2 - meanX^2;
    if varX < 0
        varX = max(varX, 0);  % numerical guard
    end
    stdX  = sqrt(varX);
    cvX   = stdX / meanX;

    % ---- Set up struct ----
    dist.alpha  = alpha;
    dist.lambda = lambda;
    dist.xmin   = xmin;

    dist.mean = meanX;
    dist.var  = varX;
    dist.cv   = cvX;

    dist.pdf = @pdf_fun;
    dist.cdf = @cdf_fun;
    dist.ppf = @ppf_fun;
    dist.rnd = @rnd_fun;

    %================= Nested functions (capture alpha, lambda, xmin, C, etc.) =============%

    function y = pdf_fun(x)
        % PDF: f(x) = C * x^(-alpha) * exp(-lambda * x), x >= xmin
        y = zeros(size(x));
        mask = (x >= xmin);
        if any(mask)
            xm = x(mask);
            y(mask) = C .* (xm .^ (-alpha)) .* exp(-lambda .* xm);
        end
    end

    function F = cdf_fun(x)
        % CDF: F(x) = 0, x < xmin
        %           1 - Γ(1-α, λx) / Γ(1-α, λxmin), x >= xmin
        F = zeros(size(x));
        mask = (x >= xmin);
        if any(mask)
            z  = lambda .* x(mask);
            Gx = upper_incomplete_gamma(s1, z);
            F(mask) = 1 - (Gx ./ G1);
        end
    end

    function x = ppf_fun(u)
        % Quantile (inverse CDF): x(u) s.t. F(x) = u.
        %
        % Solve:
        %   1 - Γ(1-α, λx) / Γ(1-α, λxmin) = u
        % => Γ(1-α, λx) = (1-u) * Γ(1-α, λxmin)
        %
        % We solve for z = λx via fzero, then x = z / lambda.

        if any(u < 0 | u > 1)
            error('u must be in [0,1].');
        end

        x = zeros(size(u));
        % edges:
        x(u <= 0) = xmin;        % F(xmin) = 0
        x(u >= 1) = xmin * 1e6;  % approximate "∞"

        mask  = (u > 0) & (u < 1);
        u_vec = u(mask);

        if isempty(u_vec)
            return;
        end

        % Sort for more efficient bracketing
        [u_sorted, idx_sorted] = sort(u_vec);
        z_sol_sorted = zeros(size(u_sorted));

        % Initial upper bracket for z
        z_lo_global = zmin;
        z_hi_global = zmin * 2;
        if z_hi_global <= zmin
            z_hi_global = zmin + 1;
        end

        for k = 1:numel(u_sorted)
            ui = u_sorted(k);

            target = (1 - ui) * G1;  % RHS: Γ(1-α, λxmin) * (1-u)

            fun = @(z) upper_incomplete_gamma(s1, z) - target;

            % Lower bracket:
            z_lo = z_lo_global;
            f_lo = fun(z_lo);

            % Find an upper bracket with opposite sign
            z_hi = z_hi_global;
            f_hi = fun(z_hi);
            iter = 0;
            while f_lo * f_hi > 0 && iter < 50
                z_hi = z_hi * 2;
                f_hi = fun(z_hi);
                iter = iter + 1;
            end

            if f_lo * f_hi > 0
                % Failed to bracket: store NaN
                z_root = NaN;
            else
                z_root = fzero(fun, [z_lo, z_hi]);
            end

            z_sol_sorted(k) = z_root;

            % Update global hi for next iterations (larger u => larger x)
            if isfinite(z_root)
                z_hi_global = max(z_hi_global, z_root * 1.5);
            end
        end

        % Put back into original order
        z_sol = zeros(size(u_vec));
        z_sol(idx_sorted) = z_sol_sorted;

        x(mask) = z_sol ./ lambda;
        x = max(x, xmin);  % numerical safety
    end

    function x = rnd_fun(n)
        % Draw n samples by inverse transform sampling
        if nargin < 1
            n = 1;
        end
        u = rand(n,1);
        x = ppf_fun(u);
    end

end

%================= Helper: upper incomplete gamma Γ(s,z) with s<=0 OK =====================%

function G = upper_incomplete_gamma(s, z)
% UPPER_INCOMPLETE_GAMMA  Compute unregularized upper incomplete gamma Γ(s,z)
% for real s (except non-positive integers).
%
% Definition:
%   Γ(s,z) = ∫_z^∞ t^(s-1) e^(-t) dt
%
% MATLAB's gammainc(z,s,'upper') expects s >= 0. For s > 0:
%
%   Q(s,z) = gammainc(z, s, 'upper') = Γ(s,z) / Γ(s)
%   => Γ(s,z) = Q(s,z) * Γ(s).
%
% For s <= 0 (non-integer), use the recurrence:
%
%   Γ(s+1, z) = s Γ(s, z) + z^s e^(-z)
%   => Γ(s, z) = (Γ(s+1, z) - z^s e^(-z)) / s
%
% We shift to s_shift = s + k > 0, compute Γ(s_shift,z), then step down.

    % Poles at s = 0, -1, -2, ...
    if abs(s - round(s)) < 1e-12 && round(s) <= 0
        G = NaN;
        return;
    end

    if s > 0
        % Safe region
        Q = gammainc(z, s, 'upper');  % regularized upper
        G = Q * gamma(s);             % unregularized Γ(s,z)
    else
        % Shift up
        k = ceil(1 - s);              % integer so that s + k >= 1
        s_shift = s + k;

        % Compute Γ(s_shift, z) with positive shape
        Q_shift = gammainc(z, s_shift, 'upper');
        G_shift = Q_shift * gamma(s_shift);

        % Step down k times:
        % Γ(s_cur, z) = (Γ(s_cur+1, z) - z^(s_cur) e^(-z)) / s_cur
        for m = 1:k
            s_cur = s_shift - 1;
            G_cur = (G_shift - z.^(s_cur) .* exp(-z)) ./ s_cur;

            G_shift = G_cur;
            s_shift = s_cur;
        end

        G = G_shift;
    end
end
