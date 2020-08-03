function N=FluxCalcsSimple(N)
%FluxCalcsSimple.m - Function to estimate surface mass balance
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
    N.Smean=N.umult.*N.S;
    N.Umean=N.umult.*N.U;
    N.Vmean=N.umult.*N.V;
    
    %% calculate flux divergence per pixel
    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));
        
    %initialize
    N.FDIV = zeros(size(N.THX));
    N.FDIVx = zeros(size(N.THX));
    N.FDIVy = zeros(size(N.THX));

    %pixel-based flux magnitude
    N.FLUX = N.Smean.*N.THX;
    N.FLUX(N.FLUX<=0)=NaN;
    
    %first order centered-difference
    N.FDIVx(:,2:end-1) = (N.Umean(:,3:end).*N.THX(:,3:end)-N.Umean(:,1:end-2).*N.THX(:,1:end-2))/2./dx;
    N.FDIVy(2:end-1,:) = (N.Vmean(3:end,:).*N.THX(3:end,:)-N.Vmean(1:end-2,:).*N.THX(1:end-2,:))/2./dy;
    N.FDIV = N.FDIVx+N.FDIVy; %total, m/yr

    %trim to mask
    N.FDIV(N.MASK==0)=NaN;

    %% aggregate variables over zones
    N.z2fdiv = zonal_aggregate(N.zones,N.FDIV); % aggregates values in the zone - simple mean excluding NaNs; same result as perimeter integration
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
   
     %% SMB
    
    N.SMB = N.Hdensity.*N.DH+N.Qdensity.*N.FDIV; %continuity equation. note that 'density' terms are actually specific gravity
    N.SMBz2= zonal_aggregate(N.zones,N.SMB); %aggregates values in the zone - simple mean

    %mask before plotting
    N.DH((N.MASK==0))=NaN;
    N.U((N.MASK==0))=NaN;
    N.V((N.MASK==0))=NaN;
    N.THX((N.MASK==0))=NaN;
    N.FDIV((N.MASK==0))=NaN;
    N.FDIVx((N.MASK==0))=NaN;
    N.FDIVy((N.MASK==0))=NaN;
    N.SMB((N.MASK==0))=NaN;
    N.zDH((N.MASK==0))=NaN;
    N.SMBz2((N.MASK==0))=NaN;
    N.z2DH((N.MASK==0))=NaN;
    N.z2fdiv((N.MASK==0))=NaN;