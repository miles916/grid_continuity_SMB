function [Bx,By] = distfilter_gradient(N,A,W_in)
%distfilter_gradient - implements a variable-range determination of x and y
%gradients using a 2D local plane fit (weighted least squares) 
%over distances according to 4x 'local' ice thickness, with optional weights.
%
% Syntax:  [Bx,By] = distfilter_gradient(N,A,W_in)
%
% Inputs:
%    N - data stack in grid_continuity_SMB structure
%    A - image for gradients with variable distance
%    W_in - (optional) weight map for pixels (0 to 1)
%
% Outputs:
%    Bx - gradient of A in x direction
%    By - gradient of A in y direction
%
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2023; Last revision: 05-Nov-2026

% Handle weights
if nargin < 3 || isempty(W_in)
    W_in = double(~isnan(A));
else
    % Ensure valid weights
    W_in(isnan(W_in)) = 0;
    W_in(N.MASK==0) = 0;
end

% construct special filter
fs.L = 4; %length of filter area in ice thicknesses
fs.sdist=ceil(fs.L.*N.THX./N.DX); %local thickness in pixel-equivalent
fs.smax=nanmax(fs.sdist(:)); %maximum distance in pixel-equivalent

% determine maximum value in radius (smoothing the distance field)
winLength=ceil(250./N.DX); %find within 250m window
fs.sdist = colfilt(fs.sdist, [winLength winLength], 'sliding', @max).*N.MASK;

% Variable ranges
STSIZE = int16(fs.sdist); STSIZE(STSIZE<=1)=2; %minimum distance for sliding filter
STSIZE(N.MASK==0) = 0;

Bx = nan(size(A));
By = nan(size(A));

% Group pixels by window size for efficiency
unique_sizes = double(unique(STSIZE(N.MASK)));

% Prepare mask and data
W = W_in; 
A_val = A;
A_val(isnan(A)) = 0;

for sz_radius = unique_sizes'
    % Window size is sz_radius*2 + 1
    sz = sz_radius * 2 + 1;
    k = sz_radius;
    
    % Relative coordinates for the window
    [x_grid, y_grid] = meshgrid(-k:k, -k:k);
    
    % Kernels for weighted sums
    % Weights W are applied to the A_val and the coordinate kernels
    
    % Compute sums using convolution with coordinate kernels
    S_w  = conv2(W, ones(sz), 'same');
    S_x  = conv2(W, x_grid, 'same');
    S_y  = conv2(W, y_grid, 'same');
    S_xx = conv2(W, x_grid.^2, 'same');
    S_yy = conv2(W, y_grid.^2, 'same');
    S_xy = conv2(W, x_grid .* y_grid, 'same');
    
    S_z  = conv2(A_val .* W, ones(sz), 'same');
    S_xz = conv2(A_val .* W, x_grid, 'same');
    S_yz = conv2(A_val .* W, y_grid, 'same');
    
    % Solve 3x3 system for each pixel in this group
    idx = find(STSIZE == sz_radius);
    for i = 1:length(idx)
        ii = idx(i);
        if S_w(ii) < 3, continue; end % Need at least 3 points for a plane
        
        M = [S_xx(ii), S_xy(ii), S_x(ii); ...
             S_xy(ii), S_yy(ii), S_y(ii); ...
             S_x(ii),  S_y(ii),  S_w(ii)];
        
        RHS = [S_xz(ii); S_yz(ii); S_z(ii)];
        
        % Check rank to avoid singular matrix
        if rank(M) == 3
            sol = M \ RHS;
            Bx(ii) = sol(1) / N.DX;
            By(ii) = -sol(2) / N.DX; % Negative because y increases downwards in image coords
        end
    end
end

Bx(N.MASK==0) = NaN;
By(N.MASK==0) = NaN;

end
