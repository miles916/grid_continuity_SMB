function N=FluxCalcsUncertainty(N)
%FluxCalcsUncertainty.m - Function to estimate surface mass balance
% distribution for a glacier from inputs of ice thickness, thinning, and
% velocity, based on the continuity equation (see, e.g. Bisset et al, 2020: https://doi.org/10.3390/rs12101563)
% In addition to inputs for FluxCalcsSimple, requires N to contain variables of the following as single value or grid in variable units:
%       sigUmult - uncertainty of column-averaged velocity factor
%       sigU - uncertainty of u component of velocity (m/a)
%       sigV - uncertainty of v component of velocity (m/a)
%       sigTHX - uncertainty of ice thickness in terms of systematic bias (percent)
%       sigTHX2 - uncertainty of ice thickness in terms of random error (m)
%       sigDH - uncertainty of thinning (m/a)
% The code then runs through nr randomly-chosen parameter sets.
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 16-June-2020

%% baseline calculations for unperturbed case
    N.Smean=N.umult.*N.S;
    N.Umean=N.umult.*N.U;
    N.Vmean=N.umult.*N.V;
    
    %pixel-based flux magnitude
    N.FLUX = N.Smean.*N.THX;N.FLUX(N.FLUX<=0)=NaN;
    
%%  nr random draws for systematic (THXscale,umult,Vscale) and random (10m thx, sig m/a vel, sig density & dH) uncertainty
%     nr= 1000;
    nr = N.uncertainty;
% Accumulate Monte Carlo moments without storing all full-grid realizations.
    sumFDIV=zeros(size(N.THX));
    sumSMB=zeros(size(N.THX));
    sumHdens=zeros(size(N.THX));
    ssFDIV=zeros(size(N.THX));
    ssSMB=zeros(size(N.THX));
    ssHdens=zeros(size(N.THX));
    countFDIV=zeros(size(N.THX));
    countSMB=zeros(size(N.THX));
    countHdens=zeros(size(N.THX));
    
%     loop through runs in parallel (parfor reductions; runs serially if no pool/toolbox)
    useProgress = license('test','Distrib_Computing_Toolbox');
    pool = [];
    if useProgress
        pool = gcp('nocreate');
    end
    nWorkers = 1;
    if ~isempty(pool)
        nWorkers = pool.NumWorkers;
    end
    if useProgress
        progressQueue = parallel.pool.DataQueue;
        afterEach(progressQueue, @(~) mcProgress(nr, nWorkers, false));
    end
    mcProgress(nr, nWorkers, true);
    tic
    parfor irun=1:nr
        C=N;

        %	determine MC adjustments
        C.umult=N.umult+randn.*N.sigUmult;
        C.THXmult=1+randn.*N.sigTHX; %systematic thickness scaling uncertainty
        C.THXrand=randn(size(C.THX)).*N.sigTHX2; %extra to determine random adjustments
        C.Uadd=randn(size(C.U)).*N.sigU;
        C.Vadd=randn(size(C.V)).*N.sigV;
        C.dHadd=randn(size(C.DH)).*N.sigDH;
        
        %   determine current inputs to flux calcs
        C.U=C.U+C.Uadd;
        C.V=C.V+C.Vadd;
        C.S=hypot(C.U,C.V);
        C.THX=(C.THX.*C.THXmult+C.THXrand).*C.MASK;
        C.DH=C.DH+C.dHadd;
        
        % calculate fluxes and SMB
        C=FluxCalcsSimple(C); 
        
        % Accumulate per-pixel sums while ignoring masked or invalid values.
        [v,x] = validValues(C.FDIV);
        countFDIV = countFDIV + v;
        sumFDIV = sumFDIV + x;
        ssFDIV = ssFDIV + x.^2;
        [v,x] = validValues(C.SMB);
        countSMB = countSMB + v;
        sumSMB = sumSMB + x;
        ssSMB = ssSMB + x.^2;
        [v,x] = validValues(C.Hdensity);
        countHdens = countHdens + v;
        sumHdens = sumHdens + x;
        ssHdens = ssHdens + x.^2;
        if useProgress
            send(progressQueue, irun);
        end
    end
    toc
    meanFDIV = sumFDIV./countFDIV;  meanFDIV(countFDIV==0) = 0;
    meanSMB = sumSMB./countSMB;  meanSMB(countSMB==0) = 0;
    meanHdens = sumHdens./countHdens;  meanHdens(countHdens==0) = 0;
    m2FDIV = ssFDIV - countFDIV.*meanFDIV.^2;
    m2SMB = ssSMB - countSMB.*meanSMB.^2;
    m2Hdens = ssHdens - countHdens.*meanHdens.^2;
    %% postprocess MC stack into N
    N.SMB=meanSMB;
    N.FDIV=meanFDIV;
    N.Hdensity=meanHdens;

    N.SMBu=sqrt(max(m2SMB,0)./max(countSMB-1,1));
    N.FDIVu=sqrt(max(m2FDIV,0)./max(countFDIV-1,1));
    N.SMBu(countSMB<2)=NaN;
    N.FDIVu(countFDIV<2)=NaN;
    N.SMB(~N.MASK)=NaN;
    N.FDIV(~N.MASK)=NaN;
    N.Hdensity(~N.MASK)=NaN;
    N.SMBu(~N.MASK)=NaN;
    N.FDIVu(~N.MASK)=NaN;
    
    N.Qdensity = 0.9; %900 kg m3 everywhere; 0.9 is actually the specific gravity

    N.z2fdiv = zonal_aggregate(N.zones,N.FDIV); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
    N.z2DH = zonal_aggregate(N.zones,N.DH); % aggregates values in the zone - simple mean
    N.SMBz2= zonal_aggregate(N.zones,N.SMB); %aggregates values in the zone - simple mean
    
    N.z2fdiv_unc = zonal_aggregate_v2(N.zones,N.FDIVu,'rssn'); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
    N.z2DH_unc = zonal_aggregate_v2(N.zones,N.sigDH,'rssn'); % aggregates values in the zone - simple mean
    N.SMBz2_unc= zonal_aggregate_v2(N.zones,N.SMBu,'rssn'); %aggregates values in the zone - simple mean
end

function mcProgress(nr, nWorkers, reset)
% Print Monte Carlo progress as persistent lines (scrollable) with an ETA.
% ETA: after the first run completes, extrapolate from its duration assuming
% nWorkers runs execute concurrently; once all workers have reported at
% least once, switch to the measured throughput.
persistent done t0 tFirst nextPct
if reset
    done = 0;
    t0 = tic;
    tFirst = NaN;
    nextPct = 10;
    fprintf('  MC: starting %d runs on %d worker(s)...\n', nr, nWorkers);
    return
end
done = done + 1;
elapsed = toc(t0);
if done == 1
    tFirst = elapsed;
    est = tFirst*ceil(nr/nWorkers);
    fprintf('  MC: first run done in %s -> est. total %s (est. remaining %s)\n', ...
        fmtSec(tFirst), fmtSec(est), fmtSec(max(est - elapsed, 0)));
    return
end
if done >= nWorkers
    remaining = (elapsed - tFirst)./(done - 1).*(nr - done);
else
    remaining = max(tFirst*ceil(nr/nWorkers) - elapsed, 0);
end
pct = 100*done/nr;
if pct >= nextPct || done == nr
    fprintf('  MC: %d/%d runs (%.0f%%) | elapsed %s | est. remaining %s\n', ...
        done, nr, pct, fmtSec(elapsed), fmtSec(remaining));
    nextPct = (floor(pct/10) + 1)*10;
end
end

function s = fmtSec(sec)
sec = round(sec);
if sec >= 3600
    s = sprintf('%dh%02dm%02ds', floor(sec/3600), floor(mod(sec,3600)/60), mod(sec,60));
elseif sec >= 60
    s = sprintf('%dm%02ds', floor(sec/60), mod(sec,60));
else
    s = sprintf('%ds', sec);
end
end

function [valid,x] = validValues(sample)
valid = isfinite(sample);
x = zeros(size(sample));
x(valid) = sample(valid);
end
