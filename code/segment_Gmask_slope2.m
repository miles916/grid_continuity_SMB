function zones = segment_Gmask_slope2(DEM,mask,dx,DL)
% DEM is DEM for whole area
% mask denotes glacier area

DEM0=DEM;

DEM1 = DEM0;
DEM1(mask==0)=NaN;

DEM2 = inpaint_nans(DEM1,0);

DEM3 = imgaussfilt(DEM2,2);
DEM3 = imgaussfilt(DEM3,2);

% mask = isnan(DEM);

DEM = DEM3;
DEM(mask==0)=NaN;

% [SLO,DIR] = imgradient(DEM./dx);
[Fx,Fy] = gradient(DEM3);
SLO = sqrt(Fx.^2+Fy.^2)/dx;
SLO(mask==0)=NaN;

% %% elevation range 
% 
ELmax = nanmax(DEM(DEM(:)>2000));
ELmin = nanmin(DEM((DEM(:)>2000)));

%% find slope vs elevation and optimal thresholds
ELs=ELmin:25:ELmax;
mSLO=0.*ELs;
for iEL=1:length(ELs)-1
    c=(DEM<ELs(iEL+1))&(DEM>=ELs(iEL));
    mSLO(iEL)=nanmedian(SLO(c));
end

% figure;
% plot(ELs(1:end-1),mSLO(1:end-1))

%% set contour intervals for each 10-centile
% ELAest = nanmedian(DEM(DEM(:)>0)); %median elevation
% ELAest = prctile(DEM(DEM(:)>0),40); %40th centile elevation
ELs = prctile(DEM(DEM(:)>0),[0:10:100]); %40th centile elevation

% DL=200;
ELvs=[];
for iEL=1:length(ELs)-1
    c=(DEM<ELs(iEL+1))&(DEM>=ELs(iEL));
    mSLO(iEL)=nanmedian(SLO(c));
    cCINT=mSLO(iEL).*DL;
    cELvs = cCINT.*[(ELs(iEL)/cCINT):1:(ELs(iEL+1)/cCINT)];
    ELvs=[ELvs,cELvs];
end

zones = zeros(size(DEM));
%% iterate through zones
iZ = 0;

for iEL=2:length(ELvs)
    cur = (DEM>ELvs(iEL-1))&(DEM<=ELvs(iEL));
    cur = imfill(cur,'holes');
    cz = bwlabel(cur); %label each segment individually
    cn = max(cz(:)); %number of new zones
    zones = zones+cz+cur.*iZ;
    iZ = iZ+cn;
end
zones=uint16(zones);


zones = remove_small_zones(zones,mask,floor(30.*(dx/50.^2)));

figure;
imagesc(zones.*uint16(mask)); colorbar;
axis square
% set(gca,'ydir','normal')
colormap([0,0,0;lines(256)])
% caxis([0,length(zonelist)])
% linkaxes([a1,a2],'xy')

