function N=FluxCalcsUncertainty(N,umult,sigUmult,sigU,sigV,sigTHX,sigDH)
%FluxCalcsUncertainty.m - Function to estimate surface mass balance
% distribution for a glacier from inputs of ice thickness, thinning, and
% velocity, based on the continuity equation (see, e.g. Bisset et al, 2020: https://doi.org/10.3390/rs12101563)
%
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 16-June-2020

%% FLUX AND EMERGENCE CALCS 
    N.Smean=umult.*N.S;
    N.Umean=umult.*N.U;
    N.Vmean=umult.*N.V;
    
    N.Uunc=sqrt((sigUmult./umult).^2+sigU.^2).*N.MASK;
    N.Vunc=sqrt((sigUmult./umult).^2+sigV.^2).*N.MASK;
    
    %% calculate flux divergence per pixel
    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));
        
    %initialize
    N.FDIV = zeros(size(N.THX));
    N.FDIVx = zeros(size(N.THX));
    N.FDIVy = zeros(size(N.THX));
    N.FDIVxU = zeros(size(N.THX));
    N.FDIVyU = zeros(size(N.THX));

    %pixel-based flux magnitude
    N.FLUX = N.Smean.*N.THX;N.FLUX(N.FLUX<=0)=NaN;
    N.FLUXx= N.Umean.*N.THX;N.FLUXy= N.Vmean.*N.THX;
    N.FLUXy= N.Umean.*N.THX;N.FLUXy= N.Vmean.*N.THX;
    N.FLUXux= sqrt(N.Uunc.^2+sigTHX.^2).*N.MASK; %normalized
    N.FLUXuy= sqrt(N.Vunc.^2+sigTHX.^2).*N.MASK; %normalized
    
    %first order centered-difference
    N.FDIVx(:,2:end-1) = (N.Umean(:,3:end).*N.THX(:,3:end)-N.Umean(:,1:end-2).*N.THX(:,1:end-2))/2./dx;
    N.FDIVy(2:end-1,:) = (N.Vmean(3:end,:).*N.THX(3:end,:)-N.Vmean(1:end-2,:).*N.THX(1:end-2,:))/2./dy;
    N.FDIV = N.FDIVx+N.FDIVy; %total, m/yr

    N.FDIVxU(:,2:end-1) = sqrt((N.FLUXux(:,3:end).*N.FLUXx(:,3:end)).^2+(N.FLUXux(:,1:end-2).*N.FLUXx(:,1:end-2)).^2)/2./dx;
    N.FDIVyU(2:end-1,:) = sqrt((N.FLUXuy(3:end,:).*N.FLUXy(3:end,:)).^2+(N.FLUXuy(1:end-2,:).*N.FLUXy(1:end-2,:)).^2)/2./dy;
    N.FDIVu = sqrt(N.FDIVxU.^2+N.FDIVyU.^2); %total, m/yr

    %trim to mask
    N.FDIV(N.MASK==0)=NaN;

    %% aggregate variables over zones
    N.z2fdiv = zonal_aggregate(N.zones,N.FDIV); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration; m/a
    N.z2fdiv_unc = zonal_aggregate_v2(N.zones,N.FDIVu,'rssn'); % aggregates values in the zone - simple rss/n excluding NaNs; same result as perimeter integration; m/a
    N.z2DH = zonal_aggregate(N.zones,N.DH); % aggregates values in the zone - simple mean

    %% density corrections
    N.Qdensity = 0.9; %900 kg m3 everywhere; 0.9 is actually the specific gravity
    N.Hdensity = 0.9.*N.MASK; %initial value everywhere of 900kg m3; 0.9 is actually the specific gravity
    
%     %zonal implementation
%     ind1=(N.z2fdiv>0); %positive emergence
%     ind2=(N.z2DH>0); %thickening
%     ind3=(abs(N.z2DH)>abs(N.z2fdiv)); %elevation change greater than fdiv
    
    %grid implementation of density correction
    ind1=(N.FDIV>0); 
    ind2=(N.DH>0);
    ind3=(abs(N.DH)>abs(N.FDIV));
    ind4 = N.DEM>nanmedian(N.DEM(N.MASK));
    
    N.Hdensity(ind1&~ind2)=0.9; %thinning and emergence = melt
    N.Hdensity(~ind1&ind2)=0.6; %thickening and submergence = acc
    N.Hdensity(ind1&ind2&ind3)=0.6; %emergence and thickening, more thickening - acc
    N.Hdensity(ind1&ind2&~ind3)=0.85; %emergence and thickening, more thickening - mixed
    N.Hdensity(~ind1&~ind2&ind3)=0.9; %submergence and thinning, more thinning - melt
    N.Hdensity(~ind1&~ind2&~ind3)=0.85; %submergence and thinning, less thinning - mixed
    N.Hdensity((ind4==0)&N.MASK)=0.9; % below median elevation
   
    N.sigDens=0.06; %uncertainty in terms of specific gravity, per Huss 2013
    %% SMB
    
    N.SMB = N.Hdensity.*N.DH+N.Qdensity.*N.FDIV; %continuity equation. note that 'density' terms are actually specific gravity
    N.SMBu = sqrt((N.Hdensity.*N.DH).^2.*(sqrt((N.sigDens./N.Hdensity).^2+(sigDH./N.DH).^2))+(N.Qdensity.*N.FDIV).^2.*(sqrt((N.sigDens./N.Qdensity).^2+(N.FDIVu./N.FDIV).^2)));
    N.SMBz2= zonal_aggregate(N.zones,N.SMB); %aggregates values in the zone - simple mean
    N.SMBz2e= zonal_aggregate_v2(N.zones,N.SMB,'rssn'); %aggregates values in the zone - simple rss/n

    %mask before plotting
    N.DH((N.MASK==0))=NaN;
    N.U((N.MASK==0))=NaN;
    N.V((N.MASK==0))=NaN;
    N.THX((N.MASK==0))=NaN;
    N.FDIV((N.MASK==0))=NaN;
    N.FDIVu((N.MASK==0))=NaN;
    N.FDIVx((N.MASK==0))=NaN;
    N.FDIVy((N.MASK==0))=NaN;
    N.SMB((N.MASK==0))=NaN;
    N.SMBu((N.MASK==0))=NaN;
    N.zDH((N.MASK==0))=NaN;
    N.SMBz2((N.MASK==0))=NaN;
    N.SMBz2e((N.MASK==0))=NaN;
    N.z2DH((N.MASK==0))=NaN;
    N.z2fdiv((N.MASK==0))=NaN;
    N.z2fdiv_unc((N.MASK==0))=NaN;