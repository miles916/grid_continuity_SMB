function N=FluxCalcsUncertainty(N,nr) 
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
    rns=randn(nr,6); %values
    prns=cdf('Normal',rns,1,1);
        
% initialize vars
    tFDIV=zeros([size(N.THX),nr]);
    tSMB=zeros([size(N.THX),nr]);
    tHdens=zeros([size(N.THX),nr]);
    
    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));
    
%     loop through runs
    tic
    for irun=1:nr
% irun =1
        C=N;

        %	determine MC adjustments
        C.umult=icdf('normal',prns(irun,1),N.umult,N.sigUmult);
        C.THXmult=icdf('normal',(prns(irun,2)),1,N.sigTHX); %extra to determine scaling up or down
        C.THXrand=randn(size(C.THX)).*N.sigTHX2; %extra to determine random adjustments
        C.Uadd=icdf('normal',(prns(irun,3)),0,N.sigU); %extra to determine random adjustments
        C.Vadd=icdf('normal',(prns(irun,4)),0,N.sigV); %extra to determine random adjustments
        C.dHadd=icdf('normal',(prns(irun,5)),0,N.sigV); %extra to determine random adjustments
        
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
        
        clear C
    end
    toc
    %% postprocess MC stack into N
    N.SMB=nanmean(tSMB,3);
    N.FDIV=nanmean(tFDIV,3);
    N.Hdensity=nanmean(tHdensity,3);
    
    N.SMBu=nanstd(tSMB,[],3);
    N.FDIVu=nanstd(tFDIV,[],3);
    N.Hdensityu=nanstd(tHdensity,[],3);
    
    N.Qdensity = 0.9; %900 kg m3 everywhere; 0.9 is actually the specific gravity

    N.z2fdiv = zonal_aggregate(N.zones,N.FDIV); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
    N.z2DH = zonal_aggregate(N.zones,N.DH); % aggregates values in the zone - simple mean
    N.SMBz2= zonal_aggregate(N.zones,N.SMB); %aggregates values in the zone - simple mean
    
    N.z2fdiv_unc = zonal_aggregate_v2(N.zones,N.FDIVu,'rssn'); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
    N.z2DH_unc = zonal_aggregate_v2(N.zones,N.sigDH,'rssn'); % aggregates values in the zone - simple mean
    N.SMBz2_unc= zonal_aggregate_v2(N.zones,N.SMBu,'rssn'); %aggregates values in the zone - simple mean
