function weights = calculate_SIA_weights(N)
%calculate_SIA_weights - computes weights based on SIA consistency.
%
% R_SIA = u_obs / u_SIA.
% Weights decrease as R_SIA deviates from 1.
%
% Author: Evan Miles, 2026

% 1. Calculate driving stress: tau = rho * g * sin(alpha) * H
rho = 900; g = 9.81;
% Calculate surface slope alpha
[Fx, Fy] = gradient(N.DEM);
alpha = sqrt(Fx.^2 + Fy.^2) / N.DX;
tau = rho * g * alpha .* N.THX;

% 2. Calculate SIA velocity: u_SIA = (2A/(n+1)) * tau^n * H
% Using standard values for temperate ice
A = 1e-16; n = 3;
u_SIA = (2*A/(n+1)) .* tau.^n .* N.THX;

% 3. Calculate Ratio R_SIA
R_SIA = N.S ./ (u_SIA + 1e-6); % Add small epsilon to avoid div by zero

% 4. Calculate Weights: Gaussian around 1
sigma = 0.5; % Flexibility parameter
weights = exp(-(R_SIA - 1).^2 ./ (2 * sigma^2));

% Mask out non-glacier pixels
weights(~N.MASK) = 0;
end
