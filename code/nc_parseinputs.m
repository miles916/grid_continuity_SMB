function N = nc_parseinputs(inputFile)
%NC_PARSEINPUTS Read ContinuIX gridded NetCDF inputs into the model struct.
% Raster fields are transposed to row=y / column=x for MATLAB image tools.
% inputFile may be a path or a cell array of paths whose variables are
% disjoint pieces of one input (e.g. base_vars + glob_* files). The first
% file defines the grid; rasters from files on other grids are resampled to it.

files = cellstr(inputFile);
for iFile = 1:numel(files)
    if ~isfile(files{iFile})
        error('ContinuIX:MissingInput', 'Input NetCDF was not found: %s', files{iFile});
    end
end

variableNames = {};
varFile = containers.Map();
for iFile = 1:numel(files)
    info = ncinfo(files{iFile});
    for v = {info.Variables.Name}
        if ~isKey(varFile, v{1})
            varFile(v{1}) = files{iFile};
            variableNames{end+1} = v{1}; %#ok<AGROW>
        end
    end
end
inputFile = files{1};
requiredVariables = {'x', 'y', 'DHDT', 'VX', 'VY', 'THK', ...
    'DEM', 'BED', 'ELEVBINS', 'ICEMASK'};
missingVariables = requiredVariables(~ismember(requiredVariables, variableNames));
if ~isempty(missingVariables)
    error('ContinuIX:MissingVariables', ...
        'Input %s is missing required variables: %s', ...
        strjoin(files, ' + '), strjoin(missingVariables, ', '));
end
x = double(ncread(inputFile, 'x'));
y = double(ncread(inputFile, 'y'));
varFile('__x') = x; % base grid; rasters on other grids are resampled to it
varFile('__y') = y;
if numel(x) < 2 || numel(y) < 2 || ...
        any(~isfinite(x(:))) || any(~isfinite(y(:))) || ...
        any(diff(x) == 0) || any(diff(y) == 0) || ...
        ~(all(diff(x) > 0) || all(diff(x) < 0)) || ...
        ~(all(diff(y) > 0) || all(diff(y) < 0))
    error('ContinuIX:InvalidCoordinates', ...
        'Expected finite, strictly monotonic projected x and y coordinates.');
end
dx = median(diff(x));
dy = median(diff(y));
if max(abs(diff(x) - dx)) > 1e-6 * abs(dx) || ...
        max(abs(diff(y) - dy)) > 1e-6 * abs(dy)
    error('ContinuIX:NonUniformGrid', ...
        'Input %s has nonuniform coordinates; this method requires a regular grid.', ...
        inputFile);
end

N = struct();
N.x3 = x(:);
N.y3 = y(:);
N.DX = abs(dx);
N.THX = readRaster(varFile, 'THK');
N.DH = readRaster(varFile, 'DHDT');
N.U = readRaster(varFile, 'VX');
N.V = readRaster(varFile, 'VY');
N.DEM = readRaster(varFile, 'DEM');
N.BED = readRaster(varFile, 'BED');
N.ELEVBINS = readRaster(varFile, 'ELEVBINS');
iceMask = readRaster(varFile, 'ICEMASK');

expectedSize = [numel(y), numel(x)];
fields = {'THX', 'DH', 'U', 'V', 'DEM', 'BED', 'ELEVBINS'};
for iField = 1:numel(fields)
    if ~isequal(size(N.(fields{iField})), expectedSize)
        error('ContinuIX:GridSizeMismatch', ...
            'Variable %s in %s does not match the x/y coordinate lengths.', ...
            fields{iField}, inputFile);
    end
end
if ~isequal(size(iceMask), expectedSize)
    error('ContinuIX:GridSizeMismatch', ...
        'ICEMASK in %s does not match the x/y coordinate lengths.', inputFile);
end

N.MASK = isfinite(N.THX) & N.THX > 0 & isfinite(iceMask) & iceMask > 0;
if ~any(N.MASK(:))
    error('ContinuIX:EmptyGlacierMask', ...
        'No glacier pixels were found in %s.', inputFile);
end
N.THX(~N.MASK) = 0;
N.DH(~N.MASK) = 0;
N.U(~N.MASK) = 0;
N.V(~N.MASK) = 0;
N.S = hypot(N.U, N.V);
N.S(~N.MASK) = 0;
N.DEM(~N.MASK) = NaN;
N.BED(~N.MASK) = NaN;

% EXP01 provides gridded velocity uncertainties for some glaciers. Missing
% uncertainty data are represented as zero perturbation, never inferred.
 [N.sigU,hasSigU] = readOptionalUncertainty(varFile, variableNames, ...
    'UNCT_VX', 'VX', expectedSize);
 [N.sigV,hasSigV] = readOptionalUncertainty(varFile, variableNames, ...
    'UNCT_VY', 'VY', expectedSize);
 [N.sigDH,hasSigDH] = readOptionalUncertainty(varFile, variableNames, ...
    'UNCT_DHDT', 'DHDT', expectedSize);
 [N.sigTHX2,hasSigTHX] = readOptionalUncertainty(varFile, variableNames, ...
    'UNCT_THK', 'THK', expectedSize);
N.sigTHX = 0; % No separate systematic THK uncertainty is supplied.
N.sigU(~N.MASK) = 0;
N.sigV(~N.MASK) = 0;
N.sigDH(~N.MASK) = 0;
N.sigTHX2(~N.MASK) = 0;
N.uncertaintySources = {};
if hasSigU, N.uncertaintySources{end+1} = 'UNCT_VX'; end
if hasSigV, N.uncertaintySources{end+1} = 'UNCT_VY'; end
if hasSigDH, N.uncertaintySources{end+1} = 'UNCT_DHDT'; end
if hasSigTHX, N.uncertaintySources{end+1} = 'UNCT_THK'; end
end

function data = readRaster(varFile, variableName)
% Read a raster; if its file uses a different x/y grid than the base grid,
% linearly resample it onto the base grid (nearest outside its extent).
srcFile = varFile(variableName);
data = double(ncread(srcFile, variableName))';
xb = varFile('__x');
yb = varFile('__y');
xs = double(ncread(srcFile, 'x'));
ys = double(ncread(srcFile, 'y'));
if isequal(xs(:), xb(:)) && isequal(ys(:), yb(:))
    return
end
if xs(end) < xs(1), xs = flipud(xs(:)); data = fliplr(data); end
if ys(end) < ys(1), ys = flipud(ys(:)); data = flipud(data); end
F = griddedInterpolant({ys(:), xs(:)}, data, 'linear', 'nearest');
data = F({yb(:), xb(:)});
end

function [uncertainty,provided] = readOptionalUncertainty(varFile, variableNames, ...
        uncertaintyName, sourceName, expectedSize)
uncertainty = zeros(expectedSize);
provided = false;
if ismember(uncertaintyName, variableNames)
    uncertainty = readRaster(varFile, uncertaintyName);
    if ~isequal(size(uncertainty), expectedSize)
        error('ContinuIX:UncertaintyGridSizeMismatch', ...
            '%s in %s does not match the input grid.', uncertaintyName, varFile(uncertaintyName));
    end
    provided = true;
elseif ismember(sourceName, variableNames)
    try
        scalarUncertainty = ncreadatt(varFile(sourceName), sourceName, 'uncertainty');
    catch
        scalarUncertainty = [];
    end
    if isnumeric(scalarUncertainty) && isscalar(scalarUncertainty) && ...
            isfinite(scalarUncertainty) && scalarUncertainty >= 0
        uncertainty(:) = double(scalarUncertainty);
            provided = true;
    elseif ischar(scalarUncertainty) || (isstring(scalarUncertainty) && isscalar(scalarUncertainty))
        numericValue = str2double(scalarUncertainty);
        if isfinite(numericValue) && numericValue >= 0
            uncertainty(:) = numericValue;
            provided = true;
        end
    end
end
uncertainty(~isfinite(uncertainty) | uncertainty < 0) = 0;
end
