function nc_writeoutputs(inputFile, outputFile, N, method, P)
%NC_WRITEOUTPUTS Write a new result NetCDF (x/y coordinates, result fields,
% and metadata only; input rasters are not copied).

inputFile = cellstr(inputFile);
inputFile = inputFile{1}; % metadata/coordinates come from the first (base) file
if ~isfile(inputFile)
    error('ContinuIX:MissingInput', 'Input NetCDF was not found: %s', inputFile);
end
outputDir = fileparts(outputFile);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if isfile(outputFile)
    delete(outputFile);
end
copyCoordinatesAndMetadata(inputFile, outputFile);

writeRaster(outputFile, 'SMB', N.SMB', 'Surface mass balance', ...
    'm i.e. yr-1', 'Surface mass balance estimated from the continuity equation.');
writeRaster(outputFile, 'FDIV', N.FDIV', 'Flux divergence', ...
    'm i.e. yr-1', 'Ice flux divergence used in the continuity equation.');
writeRaster(outputFile, 'DENSITY', 1000 .* N.Hdensity', ...
    'Density', 'kg m-3', 'Effective density applied to the DHDT term.');

if isfield(N, 'SMBu') && isfield(N, 'FDIVu')
    writeRaster(outputFile, 'UNCT_SMB', N.SMBu', ...
        'Uncertainty in surface mass balance', 'm i.e. yr-1', ...
        'Standard deviation of Monte Carlo SMB estimates.');
    writeRaster(outputFile, 'UNCT_FDIV', N.FDIVu', ...
        'Uncertainty in flux divergence', 'm i.e. yr-1', ...
        'Standard deviation of Monte Carlo flux-divergence estimates.');
    ncwriteatt(outputFile, 'UNCT_SMB', 'uncertainty_method', ...
        sprintf('Standard deviation from %d Monte Carlo runs.', P.uncertainty));
    ncwriteatt(outputFile, 'UNCT_FDIV', 'uncertainty_method', ...
        sprintf('Standard deviation from %d Monte Carlo runs.', P.uncertainty));
else
    ncwriteatt(outputFile, '/', 'uncertainty_note', ...
        ['SMB and FDIV uncertainties were not estimated. The EXP01 input ' ...
         'does not provide DHDT or THK uncertainty fields.']);
end
ncwriteatt(outputFile, '/', 'continuix_uncertainty_inputs', ...
    ['UNCT_VX/UNCT_VY and numeric uncertainty attributes where provided; ' ...
     'unknown or absent DHDT/THK uncertainty treated as zero perturbation.']);
if isfield(N, 'uncertaintySources') && ~isempty(N.uncertaintySources)
    ncwriteatt(outputFile, '/', 'continuix_uncertainty_source_variables', ...
        strjoin(N.uncertaintySources, ', '));
else
    ncwriteatt(outputFile, '/', 'continuix_uncertainty_source_variables', ...
        'None supplied; velocity-multiplier uncertainty is model-derived.');
end

ncwriteatt(outputFile, '/', 'continuix_experiment', P.experiment);
ncwriteatt(outputFile, '/', 'continuix_glacier', P.Glacier);
if isfield(P, 'group')
    ncwriteatt(outputFile, '/', 'continuix_submission_group', P.group);
end
ncwriteatt(outputFile, '/', 'continuix_method', method.name);
ncwriteatt(outputFile, '/', 'continuix_zone_mode', method.zoneMode);
ncwriteatt(outputFile, '/', 'continuix_result_mode', method.resultMode);
ncwriteatt(outputFile, '/', 'continuix_fdivfilt', method.fdivfilt);
ncwriteatt(outputFile, '/', 'continuix_velocity_multiplier_setting', P.umult);
ncwriteatt(outputFile, '/', 'continuix_grid_spacing_m', N.DX);
if isfield(P, 'processingDX')
    ncwriteatt(outputFile, '/', 'continuix_processing_grid_spacing_m', P.processingDX);
end
ncwriteatt(outputFile, '/', 'continuix_velocity_filter', P.Vfilter);
ncwriteatt(outputFile, '/', 'continuix_dhdt_filter', P.dhfilter);
ncwriteatt(outputFile, '/', 'continuix_thickness_filter', P.THXfilter);
ncwriteatt(outputFile, '/', 'continuix_segment_distance_m', P.segdist);
ncwriteatt(outputFile, '/', 'continuix_uncertainty_runs', P.uncertainty);
ncwriteatt(outputFile, '/', 'continuix_plotouts', P.plotouts);
if isfield(P, 'timing')
    ncwriteatt(outputFile, '/', 'continuix_computation_time_total_s', P.timing.total_s);
    ncwriteatt(outputFile, '/', 'continuix_computation_time_method_s', P.timing.method_s);
    ncwriteatt(outputFile, '/', 'continuix_computation_time_input_prep_s', P.timing.prep_s);
end
if isfield(P, 'compute')
    c = P.compute;
    ncwriteatt(outputFile, '/', 'continuix_compute_host', c.host);
    ncwriteatt(outputFile, '/', 'continuix_compute_processor', c.cpu);
    ncwriteatt(outputFile, '/', 'continuix_compute_physical_cores', c.cores);
    ncwriteatt(outputFile, '/', 'continuix_compute_parallel_workers', c.workers);
    ncwriteatt(outputFile, '/', 'continuix_compute_ram_GB', c.ramGB);
    ncwriteatt(outputFile, '/', 'continuix_compute_software', ['MATLAB ' c.matlab ' (' c.os ')']);
end
fprintf('Wrote %s\n', outputFile);
end

function copyCoordinatesAndMetadata(inputFile, outputFile)
% Create outputFile holding only the x/y coordinates, the spatial_ref
% variable (if any), and the input's global attributes; no input rasters.
info = ncinfo(inputFile);
names = {info.Variables.Name};
keep = intersect({'x', 'y', 'spatial_ref'}, names, 'stable');
for i = 1:numel(keep)
    v = info.Variables(strcmp(names, keep{i}));
    dimArgs = {};
    for j = 1:numel(v.Dimensions)
        dimArgs = [dimArgs, {v.Dimensions(j).Name, v.Dimensions(j).Length}]; %#ok<AGROW>
    end
    if isempty(dimArgs)
        nccreate(outputFile, v.Name, 'Datatype', v.Datatype, 'Format', 'netcdf4');
    else
        nccreate(outputFile, v.Name, 'Dimensions', dimArgs, ...
            'Datatype', v.Datatype, 'Format', 'netcdf4');
    end
    ncwrite(outputFile, v.Name, ncread(inputFile, v.Name));
    for a = v.Attributes
        if ~strcmp(a.Name, '_FillValue')
            ncwriteatt(outputFile, v.Name, a.Name, a.Value);
        end
    end
end
for a = info.Attributes
    ncwriteatt(outputFile, '/', a.Name, a.Value);
end
end

function writeRaster(filename, variableName, data, longName, units, description)
if any(~isfinite(data(:)) & ~isnan(data(:)))
    error('ContinuIX:InvalidOutputValues', ...
        'Output variable %s contains Inf values.', variableName);
end
nccreate(filename, variableName, ...
    'Dimensions', {'x', size(data, 1), 'y', size(data, 2)}, ...
    'Datatype', 'double', 'FillValue', NaN);
ncwrite(filename, variableName, data);
ncwriteatt(filename, variableName, 'long_name', longName);
ncwriteatt(filename, variableName, 'description', description);
ncwriteatt(filename, variableName, 'units', units);
ncwriteatt(filename, variableName, 'grid_mapping', 'spatial_ref');
end
