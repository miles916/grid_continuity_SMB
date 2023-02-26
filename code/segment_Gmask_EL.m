function zones = segment_Gmask_EL(DEM,mask,dx,dh)
%segment_Gmask_slope2 - Function to segment a glacier mask into
%elevation bands spaced at a pseudouniform distance based on the longitudinal surface slope
%
% Syntax:  zones = segment_Gmask_EL(DEM,mask,dx,DL)
%
% Inputs:
%    DEM - Digital Elevation Model for the entire area of analysis
%    mask - binary raster of glacier area
%    dx - cell size
%    dh - real-world approximate discretization (in same units as DEM)
%
% Outputs:
%    zones - uint16 raster corresponding to the glacier segments
%
% Other m-files required: inpaint_nans.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-June-2020

%% ------------- BEGIN CODE --------------
%crop DEM to mask
DEM1 = DEM;
DEM1(mask==0)=NaN;

%fill voids in DEM with tensioned-plate interpolation
% DEM2 = inpaint_nans(DEM1,0);

%filter DEM with low-pass gaussian 2x; this removes surface undulations to
%give something more akin to longitudinal surface gradient
DEM3 = imgaussfiltNaN(DEM1,2);
DEM3 = imgaussfiltNaN(DEM3,2);


% [SLO,DIR] = imgradient(DEM./dx);
[Fx,Fy] = gradient(DEM3);
SLO = sqrt(Fx.^2+Fy.^2)/dx;
SLO(mask==0)=NaN;

% elevation range of denoised DEM
ELmax = nanmax(DEM3(DEM3(:)>2000));
ELmin = nanmin(DEM3(DEM3(:)>2000));

% %% find relationship of slope vs elevation and optimal thresholds
% ELs=ELmin:25:ELmax; %25m segmentation
% mSLO=0.*ELs;
% for iEL=1:length(ELs)-1
%     c=(DEM3<ELs(iEL+1))&(DEM3>=ELs(iEL));
%     mSLO(iEL)=nanmedian(SLO(c)); %median slope value within each section
% end
% 
% % figure;
% % plot(ELs(1:end-1),mSLO(1:end-1))

%% set contour intervals for each decile
ELs = dh.*[floor(ELmin./dh):ceil(ELmax./dh)];

%% iterate through zones and relabel
zones = zeros(size(DEM)); %intialize
iZ = 0;

for iEL=2:length(ELs) 
    cur = (DEM3>ELs(iEL-1))&(DEM3<=ELs(iEL)); %find current semgent
    cur = imfill(cur,'holes'); %remove any holes if needed
    cz = bwlabel(cur); %label each segment individually - this gives different values for different tributaries
    cn = max(cz(:)); %number of new zones
    
    zones = zones+cz+cur.*iZ; %label current zones
    iZ = iZ+cn; %advance index
end
zones=uint16(zones); %convert to uint16

% %% remove small zones, display
% zones = remove_small_zones(zones,mask,floor(30.*(dx/50.^2))); %function to remove orphaned pixels

zones= zones.*uint16(mask);

figure;
imagesc(zones.*uint16(mask)); colorbar;
axis square
% set(gca,'ydir','normal')
colormap([0,0,0;lines(256)])
% caxis([0,length(zonelist)])
% linkaxes([a1,a2],'xy')

