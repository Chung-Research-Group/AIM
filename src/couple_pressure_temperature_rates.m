function [pressure_rate, temperature_rate] = ...
        couple_pressure_temperature_rates( ...
        base_pressure_rate, base_temperature_rate, ...
        pressure, temperature, pressure_energy_coefficient)
%COUPLE_PRESSURE_TEMPERATURE_RATES Solve the local algebraic rate coupling.
%   The dimensionless non-isothermal balances contain
%
%     dT/dtau = A - c*dP/dtau
%     dP/dtau = B + (P/T)*dT/dtau
%
%   where A and B collect all explicitly evaluated temperature and pressure
%   terms. Solving these equations simultaneously avoids evaluating the
%   pressure-work contribution with an uninitialized pressure derivative.

    validateattributes(base_pressure_rate, {'numeric'}, {'real', 'finite'});
    validateattributes(base_temperature_rate, {'numeric'}, {'real', 'finite'});
    validateattributes(pressure, {'numeric'}, {'real', 'finite', 'nonnegative'});
    validateattributes(temperature, {'numeric'}, {'real', 'finite', 'positive'});
    validateattributes(pressure_energy_coefficient, {'numeric'}, ...
        {'real', 'finite', 'nonnegative'});

    target_size = size(base_pressure_rate);
    inputs = {base_temperature_rate, pressure, temperature, ...
              pressure_energy_coefficient};
    if any(cellfun(@(value) ~isequal(size(value), target_size), inputs))
        error('AIM:BreakLab:CouplingShape', ...
              'All pressure-temperature coupling inputs must have equal size.');
    end

    pressure_to_temperature = pressure ./ temperature;
    denominator = 1 + pressure_energy_coefficient .* pressure_to_temperature;
    if any(~isfinite(denominator), 'all') || any(denominator <= 0, 'all')
        error('AIM:BreakLab:InvalidCoupling', ...
              'Pressure-temperature coupling produced an invalid denominator.');
    end

    temperature_rate = (base_temperature_rate - ...
        pressure_energy_coefficient .* base_pressure_rate) ./ denominator;
    pressure_rate = base_pressure_rate + ...
        pressure_to_temperature .* temperature_rate;

    if any(~isfinite(pressure_rate), 'all') || ...
            any(~isfinite(temperature_rate), 'all')
        error('AIM:BreakLab:NonFiniteCoupledRate', ...
              'Pressure-temperature coupling produced NaN or Inf.');
    end
end
