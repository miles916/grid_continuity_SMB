function logfile = write_runlog(action,varargin)
%write_runlog - Incrementally writes a plain-text log of a grid_continuity_SMB run for one glacier.
% Each call opens the file, appends, and closes it, so the log is complete up to
% the last finished stage even if the run later fails.
%
% Syntax:
%    logfile = write_runlog('start',outdir,P,inputs)   - create log: processor setup, settings, input files
%              write_runlog('stage',logfile,name,secs) - append the timing of a finished stage
%              write_runlog('finish',logfile,outdir,P,N,T) - append grid, flux divergence mismatch, output files
%              write_runlog('error',logfile,ME)        - append the error that stopped the run
%
% Inputs:
%    outdir - output directory (the log is written here and its files are listed at the end)
%    P - settings structure from the master script
%    inputs - cell array of input file paths
%    N - results structure after the flux calculations
%    T - structure of stage timings in seconds
%    ME - MException object

%% ------------- BEGIN CODE --------------
logfile = '';
switch action
    case 'start'
        outdir = varargin{1}; P = varargin{2}; inputs = varargin{3};
        logfile = fullfile(outdir,[P.Glacier '_runlog.txt']);
        fid = openlog(logfile,'w');
        if fid<0, return, end
        cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fprintf(fid,'grid_continuity_SMB run log - %s\n',P.Glacier);
        fprintf(fid,'Started: %s\n\n',stamp());

        fprintf(fid,'== PROCESSOR SETUP ==\n');
        fprintf(fid,'MATLAB: %s\n',version);
        fprintf(fid,'Platform: %s\n',computer);
        fprintf(fid,'OS: %s\n',strtrim(evalc('disp(system_dependent(''getos''))')));
        fprintf(fid,'CPU: %s\n',getenv('PROCESSOR_IDENTIFIER'));
        fprintf(fid,'Logical processors: %s\n',getenv('NUMBER_OF_PROCESSORS'));
        fprintf(fid,'MATLAB computational threads: %d\n',maxNumCompThreads);
        try
            [~,sysmem] = memory;
            fprintf(fid,'Physical memory: %.1f GB total, %.1f GB available\n', ...
                sysmem.PhysicalMemory.Total/2^30,sysmem.PhysicalMemory.Available/2^30);
        catch
            fprintf(fid,'Physical memory: not available on this platform\n');
        end
        if license('test','Distrib_Computing_Toolbox')
            pool = gcp('nocreate');
            if isempty(pool)
                fprintf(fid,'Parallel Computing Toolbox: available, no pool open at start\n');
            else
                fprintf(fid,'Parallel Computing Toolbox: available, pool "%s" with %d workers\n',pool.Cluster.Profile,pool.NumWorkers);
            end
        else
            fprintf(fid,'Parallel Computing Toolbox: not available (serial execution)\n');
        end
        fprintf(fid,'\n== SETTINGS ==\n');
        fn = fieldnames(P);
        for i = 1:numel(fn)
            fprintf(fid,'%s = %s\n',fn{i},fmtval(P.(fn{i})));
        end
        fprintf(fid,'\n== INPUT FILES ==\n');
        for i = 1:numel(inputs)
            d = dir(inputs{i});
            if isempty(d)
                fprintf(fid,'%s  (not found)\n',inputs{i});
            else
                fprintf(fid,'%s  (%.2f MB, modified %s)\n',inputs{i},d(1).bytes/2^20,d(1).date);
            end
        end
        fprintf(fid,'\n== PROGRESS / TIMING (s) ==\n');

    case 'stage'
        logfile = varargin{1};
        fid = openlog(logfile,'a');
        if fid<0, return, end
        cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fprintf(fid,'%-28s %10.2f   (finished %s)\n',varargin{2},varargin{3},stamp());

    case 'finish'
        logfile = varargin{1}; outdir = varargin{2}; P = varargin{3}; N = varargin{4}; T = varargin{5};
        fid = openlog(logfile,'a');
        if fid<0, return, end
        cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
        if isfield(N,'runtime_MC')
            fprintf(fid,'%-28s %10.2f   (%d runs, included in flux_calcs)\n','monte_carlo_loop',N.runtime_MC,N.uncertainty);
        end
        if isfield(T,'total')
            fprintf(fid,'%-28s %10.2f\n','total',T.total);
        end
        fprintf(fid,'\n== GRID ==\n');
        fprintf(fid,'%d rows x %d columns at %g m; glacier mask %d pixels (%.3f km^2)\n',size(N.MASK,1),size(N.MASK,2),N.DX, ...
            nnz(N.MASK),nnz(N.MASK).*N.DX.^2./1e6);

        fprintf(fid,'\n== FLUX DIVERGENCE MISMATCH ==\n');
        area = nnz(N.MASK).*N.DX.^2; %m^2
        if isfield(N,'FDIV0')
            pre = mean(N.FDIV0(N.MASK),'omitnan');
            fprintf(fid,'Mean FDIV before filtering: %.6g m/a (integrated %.6g m^3/a)\n',pre,pre.*area);
        end
        post = mean(N.FDIV(N.MASK),'omitnan');
        fprintf(fid,'Mean FDIV in output:        %.6g m/a (integrated %.6g m^3/a)\n',post,post.*area);
        if isfield(N,'delFDIV') && ~isnan(N.delFDIV)
            fprintf(fid,'Filter mismatch (mean after - before): %.6g m/a (integrated %.6g m^3/a)\n',N.delFDIV,N.delFDIV.*area);
            if isfield(N,'delFDIVu')
                fprintf(fid,'  spread over Monte Carlo runs (std): %.6g m/a\n',N.delFDIVu);
            end
        else
            fprintf(fid,'Filter mismatch: not applicable (no FDIV filtering applied)\n');
        end

        fprintf(fid,'\n== OUTPUT FILES (%s) ==\n',outdir);
        d = dir(outdir);
        d = d(~[d.isdir]);
        for i = 1:numel(d)
            if ~strcmp(d(i).name,[P.Glacier '_runlog.txt'])
                fprintf(fid,'%s  (%.2f MB)\n',d(i).name,d(i).bytes/2^20);
            end
        end
        fprintf(fid,'\nCompleted: %s\n',stamp());

    case 'error'
        logfile = varargin{1}; ME = varargin{2};
        fid = openlog(logfile,'a');
        if fid<0, return, end
        cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fprintf(fid,'\n== RUN FAILED (%s) ==\n%s\n',stamp(),ME.message);
        for i = 1:numel(ME.stack)
            fprintf(fid,'  in %s (line %d)\n',ME.stack(i).name,ME.stack(i).line);
        end

    otherwise
        error('write_runlog:action','Unknown action "%s"',action)
end
end

function fid = openlog(logfile,perm)
fid = fopen(logfile,perm);
if fid<0
    warning('write_runlog:open','Could not open %s for writing',logfile);
end
end

function s = stamp()
s = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
end

function s = fmtval(v)
if ischar(v) || isstring(v)
    s = char(v);
elseif isnumeric(v) || islogical(v)
    s = mat2str(v);
else
    s = ['<' class(v) '>'];
end
end
