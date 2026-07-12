function flux_w = WENO(flux_c, FlowDir)
%WENO Apply a second-order weighted essentially non-oscillatory scheme.
%   flux_c must contain one ghost cell at each side of the physical domain.
%   FlowDir must be either 'upwind' or 'downwind'.

    if nargin ~= 2
        error('AIM:WENO:InvalidInputCount', ...
              'WENO requires flux_c and FlowDir.');
    end
    if ~isnumeric(flux_c) || ndims(flux_c) ~= 2 || size(flux_c, 1) < 4
        error('AIM:WENO:InvalidFluxArray', ...
              'flux_c must be a numeric 2-D array with at least four rows.');
    end
    if any(~isfinite(flux_c), 'all')
        error('AIM:WENO:NonFiniteFlux', ...
              'flux_c contains NaN or Inf values.');
    end

    [n_total, n_fields] = size(flux_c);
    N = n_total - 2;
    flux_w = zeros(N + 1, n_fields, 'like', flux_c);

    % Positive, scale-aware regularization for the smoothness indicators.
    % Adding epsilon to beta, rather than to the signed flux difference,
    % prevents a zero denominator when delta_flux = -epsilon.
    local_scale = max(abs(flux_c), [], 1);
    local_scale = max(local_scale, 1);
    epsilon_weno = 1e-12 .* local_scale.^2;

    flux_w(1, :) = flux_c(1, :);
    flux_w(N + 1, :) = flux_c(N + 2, :);

    if strcmpi(FlowDir, 'upwind')
        delta0 = flux_c(3:N+1, :) - flux_c(2:N, :);
        delta1 = flux_c(3:N, :) - flux_c(2:N-1, :);

        beta0 = delta0.^2;
        beta1 = zeros(N - 1, n_fields, 'like', flux_c);
        beta1(2:end, :) = delta1.^2;
        beta1(1, :) = (2 .* (flux_c(2, :) - flux_c(1, :))).^2;

        alpha0 = (2/3) ./ (beta0 + epsilon_weno).^2;
        alpha1 = (1/3) ./ (beta1 + epsilon_weno).^2;
        weight_sum = alpha0 + alpha1;

        candidate0 = (flux_c(2:N, :) + flux_c(3:N+1, :)) ./ 2;
        candidate1 = zeros(N - 1, n_fields, 'like', flux_c);
        candidate1(2:end, :) = 1.5 .* flux_c(3:N, :) - 0.5 .* flux_c(2:N-1, :);
        candidate1(1, :) = 2 .* flux_c(2, :) - flux_c(1, :);

        flux_w(2:N, :) = (alpha0 .* candidate0 + alpha1 .* candidate1) ./ weight_sum;

    elseif strcmpi(FlowDir, 'downwind')
        delta0 = flux_c(2:N, :) - flux_c(3:N+1, :);
        delta1 = flux_c(3:N, :) - flux_c(4:N+1, :);

        beta0 = delta0.^2;
        beta1 = zeros(N - 1, n_fields, 'like', flux_c);
        beta1(1:end-1, :) = delta1.^2;
        beta1(end, :) = (2 .* (flux_c(N+1, :) - flux_c(N+2, :))).^2;

        alpha0 = (2/3) ./ (beta0 + epsilon_weno).^2;
        alpha1 = (1/3) ./ (beta1 + epsilon_weno).^2;
        weight_sum = alpha0 + alpha1;

        candidate0 = (flux_c(2:N, :) + flux_c(3:N+1, :)) ./ 2;
        candidate1 = zeros(N - 1, n_fields, 'like', flux_c);
        candidate1(1:end-1, :) = 1.5 .* flux_c(3:N, :) - 0.5 .* flux_c(4:N+1, :);
        candidate1(end, :) = 2 .* flux_c(N+1, :) - flux_c(N+2, :);

        flux_w(2:N, :) = (alpha0 .* candidate0 + alpha1 .* candidate1) ./ weight_sum;

    else
        error('AIM:WENO:InvalidFlowDirection', ...
              'FlowDir must be either ''upwind'' or ''downwind''.');
    end

    if any(~isfinite(flux_w), 'all')
        error('AIM:WENO:NonFiniteOutput', ...
              'WENO reconstruction produced NaN or Inf values.');
    end
end
