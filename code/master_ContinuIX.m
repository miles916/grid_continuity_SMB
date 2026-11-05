%master_ContinuIX_EXP01.m - Master script for ContinuIX EXP01.
%Estimate glacier surface mass balance from the continuity equation.

%% ------------- BEGIN CODE --------------
clear
close all

%% Set home directory and experiment
scriptDir = fileparts(mfilename('fullpath'));
homedir = fileparts(scriptDir);
submissionDir = fileparts(homedir);
addpath(genpath(scriptDir));

P.outputRoot = fullfile(submissionDir, 'results');
experimentRoot = 'C:\Users\miles\Dropbox\collaborations\ContinuIX\Miles_submission\Experiment_inputs\';
% experimentList = setdiff(15:20, 3); % EXP03 skipped pending further testing
experimentList = 1:20;

%% Settings
P.plotouts = 1; % Output diagnostic plots.
P.plotsHidden = true; % true: render and save plots off-screen (no windows pop up or steal focus).
P.DX = []; % Empty uses native input grid resolution.
P.maxDX = 20; % Grids finer than this are block-averaged to ~this spacing (m); coarser grids are left as is.
P.segdist = 300; % Effective linear distance between flowbands (m).
P.Vfilter = 0; % Gaussian low-pass filter on velocity.
P.dhfilter = 0; % Remove 3-sigma elevation-change outliers and smooth.
P.THXfilter = 0; % Gaussian filter on thickness.
P.umult = 0; % 0 estimates column-average velocity correction.
P.uncertainty = 30; % Monte Carlo runs using available input uncertainty fields.
P.overwrite = true; % false: skip site/method outputs that already exist; true: recompute and overwrite.
P.numWorkers = 4; % Parallel workers for Monte Carlo runs (reduce if running out of memory; 0 = serial).
P.sites = {}; % Process only these sites, e.g. {'G01','S02'}; empty = all sites.
% Python used to homogenize the raw EXP02 inputs (needs numpy, scipy, rasterio, pyshp, netCDF4, matplotlib).
P.pythonExe = 'C:\Users\miles\.copilot\session-state\2f42a3bf-cb0a-43e2-aa5a-a34fe17e8e5a\files\venv\Scripts\python.exe';
P.homogDir = fullfile(P.outputRoot, 'homogenization'); % Homogenized EXP02 inputs, diagnostic plots and summaries.

% Each method writes to <outputRoot>\<method>\<experiment>\.
% Both methods use the supplied ELEVBINS as their zonal grouping.
% methods(1) = struct('name', 'Miles2021_original', 'group', 'Miles2021_original', ...
%     'zoneMode', 'elevation_bins', 'resultMode', 'elevation_bins', ...
%     'fdivfilt', 0);
% methods(2) = struct('name', 'Miles2021_updated', 'group', 'Miles2021_updated', ...
    % 'zoneMode', 'elevation_bins', 'resultMode', 'gridded', ...
    % 'fdivfilt', 1);
methods(1) = struct('name', 'Miles2021_updated', 'group', 'Miles2021_updated', ...
    'zoneMode', 'elevation_bins', 'resultMode', 'gridded', ...
    'fdivfilt', 1); % duplicated here as first (only) in list

% Start a parallel pool for the Monte Carlo runs (falls back to serial).
% Each worker holds its own copy of the glacier grids, so fewer workers
% means less memory.
if license('test', 'Distrib_Computing_Toolbox')
    pool = gcp('nocreate');
    if P.numWorkers < 1
        if ~isempty(pool), delete(pool); end
    elseif isempty(pool) || pool.NumWorkers ~= P.numWorkers
        if ~isempty(pool), delete(pool); end
        parpool(P.numWorkers);
    end
end

P.compute = getComputeInfo(P);

% Figures are created invisible and restored to the previous default on exit.
if P.plotsHidden
    prevVisible = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    restoreVisibility = onCleanup(@() set(0, 'DefaultFigureVisible', prevVisible));
end

%% Homogenize the raw EXP02 inputs (only written if missing, or if P.overwrite)
if any(experimentList == 2)
    runEXP02Homogenization(P, experimentRoot, scriptDir);
end

%% Loop over experiments
% Pre-scan so overall progress and ETA can be reported (existing outputs are
% not counted when P.overwrite is false).
totalUnits = 0;
for expNum = experimentList
    e = sprintf('EXP%02d', expNum);
    scanJobs = experimentJobs(P, experimentRoot, expNum);
    for iJob = 1:numel(scanJobs)
        for iM = 1:numel(methods)
            if P.overwrite || ~isfile(fullfile(P.outputRoot, methods(iM).group, e, [scanJobs(iJob).stem '_' methods(iM).group '.nc']))
                totalUnits = totalUnits + 1;
            end
        end
    end
end
unitsDone = 0;
runClock = tic;
fprintf('Total: %d site/method jobs across %d experiments.\n', totalUnits, numel(experimentList));

% Complete run of every experiment for one configuration before starting the
% next, each with its own log: <outputRoot>\<group>\log_<group>.txt.
allMethods = methods;
for iRun = 1:numel(allMethods)
methods = allMethods(iRun);
P.logFile = fullfile(P.outputRoot, methods.group, ['log_' methods.group '.txt']);
writeLogHeader(P, methods);
fprintf('\n################ %s (configuration %d/%d) ################\n', methods.group, iRun, numel(allMethods));

for expNum = experimentList
P.experiment = sprintf('EXP%02d', expNum);
P.inputDir = experimentInputDir(P, experimentRoot, expNum);
fprintf('\n===== %s (%s) =====\n', P.experiment, datestr(now, 'HH:MM:SS'));

if ~exist(P.inputDir, 'dir')
    warning('ContinuIX:MissingInputDirectory', ...
        'Experiment input directory was not found, skipping: %s', P.inputDir);
    continue
end
jobs = experimentJobs(P, experimentRoot, expNum);
if isempty(jobs)
    warning('ContinuIX:NoGlacierFiles', ...
        'No usable input NetCDF files for %s in %s; skipping.', ...
        P.experiment, P.inputDir);
    continue
end

%% Process each site and method
for iGlacier = 1:numel(jobs)
    inputFile = jobs(iGlacier).files; % cell array; split inputs are merged by nc_parseinputs
    inputStem = jobs(iGlacier).stem;
    P.Glacier = jobs(iGlacier).site;

    % Skip sites whose outputs all exist already (unless overwriting)
    needed = true(1, numel(methods));
    if ~P.overwrite
        for iM = 1:numel(methods)
            needed(iM) = ~isfile(fullfile(P.outputRoot, methods(iM).group, ...
                P.experiment, [inputStem '_' methods(iM).group '.nc']));
        end
        if ~any(needed)
            fprintf('\nSkipping %s: outputs already exist (P.overwrite = false)\n', inputStem);
            continue
        end
    end
    fprintf('\nProcessing %s (site %d/%d, %d input file(s))\n', ...
        inputStem, iGlacier, numel(jobs), numel(inputFile));

    glacierClock = tic;
    %% Load all input data and prepare the native common grid
    N = nc_parseinputs(inputFile);
    if ~isempty(P.DX) && abs(P.DX - N.DX) > eps(max(P.DX, N.DX))
        error('ContinuIX:GridResolutionMismatch', ...
            'Input %s does not match configured processing resolution %.3g m.', ...
            inputStem, P.DX);
    end

    %% Coarsen overly fine grids (results are returned to the native grid on export)
    [N, coarsenInfo] = coarsen_inputs(N, P.maxDX);
    if coarsenInfo.factor > 1
        fprintf('  Coarsened %g m -> %g m grid (factor %d): %d x %d -> %d x %d px\n', ...
            coarsenInfo.nativeDX, N.DX, coarsenInfo.factor, ...
            coarsenInfo.nativeSize(2), coarsenInfo.nativeSize(1), ...
            size(N.MASK, 2), size(N.MASK, 1));
    end

    %% Cleaning raw data if needed
    V.iERR = (abs(N.U) > 400) | (abs(N.V) > 400);
    N.U(V.iERR) = NaN;
    N.V(V.iERR) = NaN;
    if P.Vfilter == 1
        N.U = imgaussfilt(N.U, 5);
        N.V = imgaussfilt(N.V, 5);
    end

    if P.THXfilter == 1
        N.THX = imgaussfilt(N.THX, 4);
        N.THX(~N.MASK) = 0;
    end

    if P.dhfilter == 1
        DH.data2 = N.DH;
        glacierDH = DH.data2(N.MASK);
        threshold = 3 .* std(glacierDH, 'omitnan');
        DH.data2(abs(DH.data2) > threshold) = NaN;
        threshold = 3 .* std(DH.data2(N.MASK), 'omitnan');
        DH.data2(abs(DH.data2) > threshold) = NaN;
        N.DH = imgaussfilt(DH.data2, 1);
    end

    %% Set fields for the legacy continuity-calculation functions
    N.fdivfilt = 0;
    N.uncertainty = P.uncertainty;
    N.DH0 = N.DH;
    fillDomain = N.MASK & isfinite(N.DEM);
    if any(fillDomain(:))
        N.DH(fillDomain) = index_nanfill( ...
            N.DH(fillDomain), N.DEM(fillDomain));
    end
    unfilledGaps = N.MASK & ~isfinite(N.DH);
    if any(unfilledGaps(:))
        validDH = N.DH(N.MASK & isfinite(N.DH));
        if isempty(validDH)
            error('ContinuIX:NoDHDTForGapFill', ...
                'No finite DHDT values remain for %s.', P.Glacier);
        end
        N.DH(unfilledGaps) = median(validDH, 'omitnan');
        warning('ContinuIX:DHDTMedianGapFill', ...
            'Filled %d unresolved DHDT pixels for %s with the glacier median.', ...
            nnz(unfilledGaps), P.Glacier);
    end
    N.DH(~N.MASK) = 0;

    %% Use the supplied ELEVBINS as the zonal grouping for every method
    N.zones = elevationBinZones(N.ELEVBINS, N.MASK);
    Nbase = N;
    prepTime = toc(glacierClock); % input loading, coarsening and gap filling shared by both methods
    for iMethod = 1:numel(methods)
        N = Nbase;
        method = methods(iMethod);
        if ~needed(iMethod)
            fprintf('  Skipping %s / %s: output exists\n', P.Glacier, method.name);
            continue
        end
        fprintf('[%s | %s %d/%d | %s %d/%d] Running Monte Carlo / flux calcs...\n', ...
            P.experiment, P.Glacier, iGlacier, numel(jobs), ...
            method.name, iRun, numel(allMethods));

        methodClock = tic;
        %% Determine column-average velocity correction
        N.fdivfilt = method.fdivfilt;
        N = derive_umult(N, struct('umult', P.umult));

        %% SMB calculations
        if N.uncertainty > 1
            N = FluxCalcsUncertainty(N);
        else
            N = FluxCalcsSimple(N);
        end

        %% Apply the selected output spatial discretization
        switch method.resultMode
            case 'elevation_bins'
                if isfield(N, 'SMBu') && isfield(N, 'FDIVu')
                    N.SMBu = zonal_aggregate_v2(N.zones, N.SMBu, 'rssn');
                    N.FDIVu = zonal_aggregate_v2(N.zones, N.FDIVu, 'rssn');
                end
                N.SMB = N.SMBz2;
                N.FDIV = N.z2fdiv;
                N.Hdensity = zonal_aggregate(N.zones, N.Hdensity);
                N.Hdensity(~N.MASK) = NaN;
            case 'gridded'
                % Keep pixel-level SMB, FDIV, and uncertainty estimates.
            otherwise
                error('ContinuIX:UnknownResultMode', ...
                    'Unsupported output mode: %s', method.resultMode);
        end

        Pmethod = P;
        Pmethod.timing = struct('method_s', toc(methodClock), 'prep_s', prepTime, ...
            'total_s', toc(methodClock) + prepTime);
        %% Export one NetCDF per glacier and method
        outputDir = fullfile(P.outputRoot, method.group, P.experiment);
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end
        outputFile = fullfile(outputDir, [inputStem '_' method.group '.nc']); % EXP##_G##_<group folder>.nc
        Pmethod.group = method.group;
        Pmethod.processingDX = N.DX;
        nc_writeoutputs(inputFile, outputFile, upsample_results(N, coarsenInfo), ...
            method, Pmethod);

        %% Output diagnostic plots
        if P.plotouts == 1
            N.resultMode = method.resultMode;
            plotTitle = sprintf('%s_%gmgrid', P.Glacier, N.DX);
            previousDir = pwd;
            restoreDir = onCleanup(@() cd(previousDir));
            cd(outputDir);
            SMBplots(N, plotTitle, [P.Glacier '_' method.name]);
            clear restoreDir
            close all
        end

        appendLog(P, P.experiment, P.Glacier, method.name, Pmethod.timing, outputFile);
        unitsDone = unitsDone + 1;
        elapsed = toc(runClock);
        eta = elapsed / unitsDone * (totalUnits - unitsDone);
        fprintf(['Finished %s / %s -> %s\n' ...
            '   OVERALL %d/%d jobs (%.1f%%) | elapsed %s | est. remaining %s\n'], ...
            P.Glacier, method.name, outputFile, unitsDone, totalUnits, ...
            100*unitsDone/totalUnits, fmtDur(elapsed), fmtDur(eta));
    end
end % glaciers
fprintf('Finished %s\n', P.experiment);
end % experiments
end % configurations

if P.plotsHidden
    set(0, 'DefaultFigureVisible', prevVisible);
end
disp('ContinuIX processing finished for all experiments.');

function inputDir = experimentInputDir(P, experimentRoot, expNum)
% EXP02 is read from its homogenized files; all other experiments from the supplied inputs.
e = sprintf('EXP%02d', expNum);
if expNum == 2
    inputDir = fullfile(P.homogDir, e);
else
    inputDir = fullfile(experimentRoot, e, e);
end
end

function runEXP02Homogenization(P, experimentRoot, scriptDir)
% Calls homogenize_EXP02.py, which writes EXP02_<site>.nc (plus a diagnostic
% plot and a processing summary) to <results>\homogenization\EXP02, skipping
% existing files unless P.overwrite is true.
outDir = fullfile(P.homogDir, 'EXP02');
cmd = sprintf('"%s" -u "%s" --input-root "%s" --out-dir "%s"', P.pythonExe, ...
    fullfile(scriptDir, 'homogenize_EXP02.py'), ...
    regexprep(experimentRoot, '[\\/]+$', ''), outDir);
if P.overwrite
    cmd = [cmd ' --overwrite'];
end
if ~isempty(P.sites)
    cmd = [cmd ' --sites ' strjoin(P.sites, ' ')];
end
fprintf('\nHomogenizing EXP02 inputs (%s)\n', datestr(now, 'HH:MM:SS'));
status = system(cmd, '-echo');
if status ~= 0
    warning('ContinuIX:HomogenizationFailed', ...
        'EXP02 homogenization reported errors (exit code %d); affected sites will be skipped.', status);
end
end

function jobs = filterSites(jobs, sites)
if ~isempty(sites) && ~isempty(jobs)
    jobs = jobs(ismember({jobs.site}, sites));
end
end

function jobs = findInputJobs(inputDir, experiment)
% One job per site (G##/S##). Single-file sites: <exp>_<site>.nc.
% Split sites: <exp>_<site>_base_vars.nc plus <exp>_<site>_glob_*.nc, which
% supply disjoint variables and are merged into a single input.
jobs = struct('site', {}, 'stem', {}, 'files', {});
if ~exist(inputDir, 'dir')
    return
end
listing = dir(fullfile(inputDir, [experiment '_*.nc']));
names = sort({listing.name});
single = regexp(names, ['^' experiment '_([GS]\d+)\.nc$'], 'tokens', 'once');
for i = find(~cellfun(@isempty, single))
    site = single{i}{1};
    jobs(end+1) = struct('site', site, 'stem', [experiment '_' site], ...
        'files', {{fullfile(inputDir, names{i})}}); %#ok<AGROW>
end
base = regexp(names, ['^' experiment '_([GS]\d+)_base_vars\.nc$'], 'tokens', 'once');
for i = find(~cellfun(@isempty, base))
    site = base{i}{1};
    isPart = startsWith(names, [experiment '_' site '_']) & endsWith(names, '.nc');
    jobs(end+1) = struct('site', site, 'stem', [experiment '_' site], ...
        'files', {fullfile(inputDir, names(isPart))}); %#ok<AGROW>
end
end

function jobs = experimentJobs(P, experimentRoot, expNum)
% EXP02 (G01-G06 and S02) is read from the homogenized files; all other
% experiments from the files in <root>\EXPxx\EXPxx.
e = sprintf('EXP%02d', expNum);
jobs = filterSites(findInputJobs(experimentInputDir(P, experimentRoot, expNum), e), P.sites);
end

function s = fmtDur(sec)
sec = round(sec);
s = sprintf('%dh%02dm%02ds', floor(sec/3600), floor(mod(sec,3600)/60), mod(sec,60));
end

function zones = elevationBinZones(elevationBins, mask)
validBins = mask & isfinite(elevationBins);
if ~any(validBins(:))
    error('ContinuIX:MissingElevationBins', ...
        'No finite ELEVBINS values are available inside the glacier mask.');
end

% Snap transition/noisy values to the nearest 50 m elevation interval;
% integer labels are only an internal encoding for zonal_aggregate.
roundedBins = 50 .* round(elevationBins(validBins) ./ 50);
[~, ~, labels] = unique(roundedBins);
zones = zeros(size(mask), 'uint16');
zones(validBins) = uint16(labels);

% Some datasets have missing ELEVBINS inside the thickness-derived mask.
% Assign only those unlabelled pixels to the nearest supplied elevation bin.
missingBins = mask & ~validBins;
if any(missingBins(:))
    [~, nearestValid] = bwdist(validBins);
    zones(missingBins) = zones(nearestValid(missingBins));
    warning('ContinuIX:MissingElevationBinsFilled', ...
        ['Assigned %d glacier pixels with missing ELEVBINS to their ' ...
         'nearest valid elevation zone.'], nnz(missingBins));
end
end

function info = getComputeInfo(P)
% Processor, core count, memory and worker information for the log and .nc metadata.
info.host = getenv('COMPUTERNAME');
info.cpu = strtrim(getenv('PROCESSOR_IDENTIFIER'));
try
    [st, out] = system('powershell -NoProfile -Command "(Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)"');
    if st == 0 && ~isempty(strtrim(out)), info.cpu = strtrim(out); end
catch
end
info.cores = feature('numcores');
info.workers = max(P.numWorkers, 0);
info.ramGB = NaN;
try
    [st, out] = system('powershell -NoProfile -Command "[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB)"');
    if st == 0, info.ramGB = str2double(strtrim(out)); end
catch
end
info.matlab = version;
info.os = computer;
end

function writeLogHeader(P, method)
logDir = fileparts(P.logFile);
if ~exist(logDir, 'dir'), mkdir(logDir); end
new = ~isfile(P.logFile);
fid = fopen(P.logFile, 'a');
c = P.compute;
if new
    fprintf(fid, 'ContinuIX computation log - %s (%s)\n', method.group, method.name);
    fprintf(fid, 'Times are wall-clock seconds. prep = input loading, coarsening and gap filling; method = Monte Carlo + flux calculations; total = prep + method.\n');
end
fprintf(fid, '\n--- Session started %s ---\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Host: %s | CPU: %s | physical cores: %d | RAM: %g GB | parallel workers: %d | MATLAB %s | %s\n', ...
    c.host, c.cpu, c.cores, c.ramGB, c.workers, c.matlab, c.os);
fprintf(fid, 'Settings: Monte Carlo runs = %d, maxDX = %g m, overwrite = %d\n', P.uncertainty, P.maxDX, P.overwrite);
fprintf(fid, '%-6s %-5s %-9s %10s %10s %10s  %s\n', 'EXP', 'site', 'method', 'prep_s', 'method_s', 'total_s', 'file');
fclose(fid);
end

function appendLog(P, experiment, site, method, timing, outputFile)
fid = fopen(P.logFile, 'a');
fprintf(fid, '%-6s %-5s %-9s %10.1f %10.1f %10.1f  %s\n', experiment, site, method, ...
    timing.prep_s, timing.method_s, timing.total_s, outputFile);
fclose(fid);
end