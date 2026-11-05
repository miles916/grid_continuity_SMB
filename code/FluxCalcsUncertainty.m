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
%       GRF_Lcorr_V_DH - correlation length for U, V, DH noise (m)
%       GRF_Lcorr_THX_mult - correlation length for THX noise (multiple of mean THX)
% The code then runs through nr randomly-chosen parameter sets.
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 05-Nov-2026

%% baseline calculations for unperturbed case
    N.Smean=N.umult.*N.S;
    N.Umean=N.umult.*N.U;
    N.Vmean=N.umult.*N.V;
    
    %pixel-based flux magnitude
    N.FLUX = N.Smean.*N.THX;N.FLUX(N.FLUX<=0)=NaN;
    
%%  nr random draws for systematic (THXscale,umult,Vscale) and random (10m thx, sig m/a vel, sig density & dH) uncertainty
%     nr= 1000;
    nr = N.uncertainty;
    rns=randn(nr,6); %values
    prns=cdf('Normal',rns,1,1);
        
% initialize vars
    tFDIV=zeros([size(N.THX),nr]);
    tSMB=zeros([size(N.THX),nr]);
    tHdensity=zeros([size(N.THX),nr]);
    tFDIV0=zeros([size(N.THX),nr]);
    tdelFDIV=nan(nr,1); %mean change in FDIV from the filter, per run
    tIntCorr=nan(nr,1); %integral correction per run
    
    % Correlation lengths in pixels
    sig_pix_V_DH = N.GRF_Lcorr_V_DH / N.DX;
    sig_pix_THX = (N.GRF_Lcorr_THX_mult * nanmean(N.THX(N.MASK))) / N.DX;

    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));
    
%     loop through runs (independent, so run in parallel; runs serially if no pool/toolbox is available)
    tMC=tic;
    parfor irun=1:nr
        C=N;

        %	determine MC adjustments (systematic/scalar)
        C.umult=icdf('normal',prns(irun,1),N.umult,N.sigUmult);
        C.THXmult=icdf('normal',(prns(irun,2)),1,N.sigTHX); %extra to determine scaling up or down
        
        %   determine MC adjustments (spatially correlated or independent random fields)
        if sig_pix_V_DH > 0
            % Spatially correlated noise (GRF)
            Uadd_f = imgaussfilt(randn(size(N.THX)), sig_pix_V_DH);
            C.Uadd = Uadd_f .* (N.sigU / std(Uadd_f(N.MASK)));
            
            Vadd_f = imgaussfilt(randn(size(N.THX)), sig_pix_V_DH);
            C.Vadd = Vadd_f .* (N.sigV / std(Vadd_f(N.MASK)));
            
            dHadd_f = imgaussfilt(randn(size(N.THX)), sig_pix_V_DH);
            C.dHadd = dHadd_f .* (N.sigDH / std(dHadd_f(N.MASK)));
        else
            % Independent random noise
            C.Uadd=randn(size(N.THX)).*N.sigU;
            C.Vadd=randn(size(N.THX)).*N.sigV;
            C.dHadd=randn(size(N.THX)).*N.sigDH;
        end
        
        if sig_pix_THX > 0
            THXrand_f = imgaussfilt(randn(size(N.THX)), sig_pix_THX);
            C.THXrand = THXrand_f .* (N.sigTHX2 / std(THXrand_f(N.MASK)));
        else
            C.THXrand=randn(size(N.THX)).*N.sigTHX2;
        end
        
        %   determine current inputs to flux calcs
        C.Umean=(C.U+C.Uadd).*C.umult;
        C.Vmean=(C.V+C.Vadd).*C.umult;
        C.THX=(C.THX.*C.THXmult+C.THXrand).*C.MASK;
        C.DH=C.DH+C.dHadd;
        
        % calculate fluxes and SMB
        C=FluxCalcsSimple(C); 
        
        % index outputs into stack
        tFDIV(:,:,irun)=C.FDIV;
        tSMB(:,:,irun)=C.SMB;
        tHdensity(:,:,irun)=C.Hdensity;
        if isfield(C,'FDIV0')
            tFDIV0(:,:,irun)=C.FDIV0;
        end
        if isfield(C,'delFDIV')
            tdelFDIV(irun)=C.delFDIV;
        end
        if isfield(C,'integral_correction')
            tIntCorr(irun)=C.integral_correction;
        end
    end
    toc
    N.runtime_MC=toc(tMC);
    %% postprocess MC stack into N
    N.SMB=nanmean(tSMB,3);
    N.FDIV=nanmean(tFDIV,3);
    N.Hdensity=nanmean(tHdensity,3);
    N.FDIV0=nanmean(tFDIV0,3);
    
    N.SMBu=nanstd(tSMB,[],3);
    N.FDIVu=nanstd(tFDIV,[],3);
    N.Hdensityu=nanstd(tHdensity,[],3);
    N.delFDIV=mean(tdelFDIV,'omitnan'); %filter mismatch of the mean FDIV, mean and spread over runs
    N.delFDIVu=std(tdelFDIV,'omitnan');
    N.integral_correction=mean(tIntCorr,'omitnan');
    
    N.Qdensity = 0.9; %900 kg m3 everywhere; 0.9 is actually the specific gravity

    N.z2fdiv = zonal_aggregate(N.zones,N.FDIV); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
    N.z2DH = zonal_aggregate(N.zones,N.DH); % aggregates values in the zone - simple mean
    N.SMBz2= zonal_aggregate(N.zones,N.SMB); %aggregates values in the zone - simple mean
    
    % Ensure uncertainties are grids for zonal aggregation
    if isscalar(N.sigDH), sigDH_grid = N.sigDH .* N.MASK; else, sigDH_grid = N.sigDH; end
    
    N.z2fdiv_unc = zonal_aggregate_v2(N.zones,N.FDIVu,'rssn'); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
    N.z2DH_unc = zonal_aggregate_v2(N.zones,sigDH_grid,'rssn'); % aggregates values in the zone - simple mean
    N.SMBz2_unc= zonal_aggregate_v2(N.zones,N.SMBu,'rssn'); %aggregates values in the zone - simple mean
