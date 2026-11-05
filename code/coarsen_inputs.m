function [N, info] = coarsen_inputs(N, targetDX)
%COARSEN_INPUTS Block-average input rasters to a coarser grid spacing.
% Only coarsens (never refines). Spacing is rounded to an integer multiple of
% the native spacing. Glacier-mask-aware: fields are averaged over glacier
% pixels only, and a coarse cell is glacier if any of its native pixels is.
% info holds what is needed by upsample_results to return to the native grid.

info = struct('factor', 1, 'nativeDX', N.DX, 'nativeSize', size(N.MASK), ...
    'nativeMask', N.MASK);
f = round(targetDX / N.DX);
if f <= 1
    return
end
info.factor = f;

[ny, nx] = size(N.MASK);
cny = ceil(ny/f);
cnx = ceil(nx/f);
count = blockSum(double(N.MASK), f, cny, cnx);
valid = count > 0;

maskedFields = {'THX', 'DH', 'U', 'V', 'DEM', 'BED', 'sigU', 'sigV', ...
    'sigDH', 'sigTHX2'};
for i = 1:numel(maskedFields)
    A = N.(maskedFields{i});
    A(~N.MASK | ~isfinite(A)) = 0;
    okCount = blockSum(double(N.MASK & isfinite(N.(maskedFields{i}))), f, cny, cnx);
    B = blockSum(A, f, cny, cnx) ./ max(okCount, 1);
    if any(strcmp(maskedFields{i}, {'DEM', 'BED'}))
        B(okCount == 0) = NaN;
    else
        B(okCount == 0) = 0;
    end
    N.(maskedFields{i}) = B;
end

% ELEVBINS: mean over finite values (zones are re-snapped to 50 m later)
E = N.ELEVBINS;
okE = isfinite(E);
E(~okE) = 0;
cntE = blockSum(double(okE), f, cny, cnx);
N.ELEVBINS = blockSum(E, f, cny, cnx) ./ max(cntE, 1);
N.ELEVBINS(cntE == 0) = NaN;

N.MASK = valid;
N.THX(~valid) = 0;
N.S = hypot(N.U, N.V);
N.S(~valid) = 0;

% Coarse cell-centre coordinates
N.x3 = blockMeanVec(N.x3, f, cnx);
N.y3 = blockMeanVec(N.y3, f, cny);
N.DX = N.DX * f;
end

function S = blockSum(A, f, cny, cnx)
[ny, nx] = size(A);
P = zeros(cny*f, cnx*f);
P(1:ny, 1:nx) = A;
S = squeeze(sum(sum(reshape(P, f, cny, f, cnx), 1), 3));
S = reshape(S, cny, cnx);
end

function c = blockMeanVec(v, f, cn)
v = v(:);
n = numel(v);
P = nan(cn*f, 1);
P(1:n) = v;
c = mean(reshape(P, f, cn), 1, 'omitnan')';
end
