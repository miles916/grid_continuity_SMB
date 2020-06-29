%master_FluxCalcsSimple.m - Master script to estimate surface mass balance
%distribution for a glacier from inputs of ice thickness, thinning, and
%velocity, based on the continuity equation (see, e.g. Bisset et al, 2020: https://doi.org/10.3390/rs12101563)
%
% New glaciers require a case to be coded with paths for the full set of
% inputs. Otherwise, simply configure the first section in terms of paths
% and desired outputs
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-June-2020

%% ------------- BEGIN CODE --------------
clear 
close all

%set home directory
homedir = 'C:\Users\miles\Documents\GitHub\Flux_thickness_iteration';

%configure inputs and outputs
plotouts=1; %output plots or not
exports=1; %save geotiffs or not
DX = 200; %resolution to run calculations at 

%choose from test glaciers
Glacier = 'Rongbuk'
% Glacier = 'Matanuska'
% Glacier = 'Aletsch'

%title of output files
datatitle = ['test_' Glacier];

addpath(genpath([homedir '\code'])) %add path to related scripts

outdir = [homedir '\results\' datatitle '_' date]
mkdir(outdir)

%% configure inputs for test glaciers
switch Glacier
    case 'Rongbuk'
        V.pathx = fullfile(homedir,'data',Glacier,'HMA_G0120_vx.tif');
        V.pathy = fullfile(homedir,'data',Glacier,'HMA_G0120_vy.tif');
        DH.path = fullfile(homedir,'data',Glacier,'15.09991_dH.tif');
        THX.path = fullfile(homedir,'data',Glacier,'thickness_RGI60-15.09991_HF2012.tif');
        DEM.path = fullfile(homedir,'data',Glacier,'15.09991_AW3D.tif');
        segdist = 300; %effective linear distance between flowbands
        V.mult=1; %scale to convert input units to m/a
        V.filter=0; %switch to smooth velocity data or not - this is the sigma for imgaussfilt  (5 is a good start, corresponds to 11x11)
        dhfilter=0; %switch to peform 2x 3-sigma outlier removal (from overall dataset - only if erroneous pixels are common)
    case 'Matanuska'
        V.pathx = fullfile(homedir,'data',Glacier,'Vx_Matanuska_ITSLIVE.tif');
        V.pathy = fullfile(homedir,'data',Glacier,'Vy_Matanuska_ITSLIVE.tif');
        DH.path = fullfile(homedir,'data',Glacier,'dhdt_Matanuska_2007_2016.tif');
%         THX.path = fullfile(homedir,'data',Glacier,'RGI60-01.10557_thickness_composite.tif');
        THX.path = fullfile(homedir,'data',Glacier,'thickness_RGI60-01.10557_HF2012.tif');
        DEM.path = fullfile(homedir,'data',Glacier,'GDEM_Matanuska.tif');
        segdist=2000; %effective linear distance between flowbands
        V.mult=1; %scale to convert input units to m/a
        V.filter=0; %switch to smooth velocity data or not - this is the sigma for imgaussfilt  (5 is a good start, corresponds to 11x11)
        dhfilter=1; %switch to peform 2x 3-sigma outlier removal (from overall dataset - only if erroneous pixels are common)
    case 'Aletsch'
        V.pathx = fullfile(homedir,'data',Glacier,'v_mean-x_(md-1)_winter_2011-2017.tiff');
        V.pathy = fullfile(homedir,'data',Glacier,'v_mean-y_(md-1)_winter_2011-2017.tiff');
        DH.path = fullfile(homedir,'data',Glacier,'Aletsch_trim_Trend(2011-2019)_m-per-year_ladfit_10m-geo.tiff');
        THX.path = fullfile(homedir,'data',Glacier,'RGI60-11.01450_thickness.tif');
        DEM.path = fullfile(homedir,'data',Glacier,'DEM_aletsch_10m_utm32.tif');
        segdist=500; %effective linear distance between flowbands
        V.mult=365.24; %scale to convert input units to m/a
        V.filter=30; %switch to smooth velocity data or not - this is the sigma for imgaussfilt (5 is a good start, corresponds to 11x11)
        dhfilter=0; %switch to peform 2x 3-sigma outlier removal (from overall dataset - only if erroneous pixels are common)
end


    outtitle1 = [num2str(DX) 'mgrid'];
    cd(outdir)
    
    %% Load all input data, resample to common grid, etc
    % load THX, determine LL bounding box
    THX.data=geotiffread(THX.path); %read data
    THX.info=geotiffinfo(THX.path); %read metadata
    THX.R = THX.info.RefMatrix; %georeferencing matrix
    
        %setup grid of coordinates for thickness data
    THX.xp=1:THX.info.Width;
    THX.yp=1:THX.info.Height;
    [THX.xm,~] = pix2map(THX.R,ones(size(THX.xp)),THX.xp);
    [~,THX.ym] = pix2map(THX.R,THX.yp,ones(size(THX.yp)));
    [THX.xmG,THX.ymG] = meshgrid(THX.xm,THX.ym);
    [THX.LatG,THX.LonG] = projinv(THX.info,THX.xmG(:),THX.ymG(:));
    
        % set initial mask from thickness data
    MASK0 = THX.data>0;
    
        % filter thickness data
    THX.dataR=THX.data;
    THX.data=imgaussfilt(THX.dataR,4); %gaussian low-pass filter. Important for thickness maps derived from field data
    THX.data(MASK0==0)=0;
    
    %project glacier's bounding box into geographic coordinates (key to subset from larger datasets)
    BBoxUTM = THX.info.BoundingBox;
    BBoxLL = [min(THX.LonG(MASK0)) min(THX.LatG(MASK0));max(THX.LonG(MASK0)) max(THX.LatG(MASK0))];
    
    % load, subset velocity
    V.info = geotiffinfo(V.pathx); %easiest way to get correct projection details
    V.R = V.info.RefMatrix; %georeferencing matrix
    
        %set up velocity coordinates
    V.xp=1:V.info.Width;
    V.yp=1:V.info.Height;
    [V.xm,~] = pix2map(V.R,ones(size(V.xp)),V.xp);
    [~,V.ym] = pix2map(V.R,V.yp,ones(size(V.yp)));
    [V.PixelRegion,V.LatG,V.LonG,V.xmG,V.ymG] = subset_geo(V.info,BBoxLL);

        %read velocity data and scale and filter if needed
    V.Uraw=V.mult.*imread(V.pathx,'PixelRegion',V.PixelRegion);
    V.Vraw=V.mult.*imread(V.pathy,'PixelRegion',V.PixelRegion);
    if V.filter>0 %gaussian low-pass filter removing extreme variations
        V.Uraw=imgaussfilt(V.Uraw,V.filter);
        V.Vraw=imgaussfilt(V.Vraw,V.filter);
    end
    
        %calculate velocity end-points (reproject full velocity vector)
    V.xm2=V.xmG+double(V.Uraw); %velocity vector end-point
    V.ym2=V.ymG+double(V.Vraw); %velocity vector end-point
    if isempty(V.info.PCS)==0 %if velocity product is projected, reproject the vectors. if not, we assume it to be oriented north and with units of m/a already
        [V.Lat2,V.Lon2]=projinv(V.info,V.xm2(:),V.ym2(:)); %unprojected velocity vector end-point
    end
    
    %identify likely errors
    iERR=(abs(V.Uraw)>400)|(abs(V.Vraw)>400);
%     histogram(V.Vraw(iERR==0),100)
    
    % load, subset DEM
    try
        DEM.info=geotiffinfo(DEM.path); %read metadata
    catch
        DEM.info=imfinfo(DEM.path);
        DEM.info=DEM.info(1);
        DEM.info.RefMatrix = [0,-DEM.info.ModelPixelScaleTag(1);DEM.info.ModelPixelScaleTag(2),0;DEM.info.ModelTiepointTag(4),DEM.info.ModelTiepointTag(5)]; %create refmatrix if needed
    end
    DEM.R = DEM.info.RefMatrix;
    [DEM.PixelRegion,DEM.LatG,DEM.LonG] = subset_geo(DEM.info,BBoxLL); %subset DEM domain
    DEM.data=imread(DEM.path,'PixelRegion',DEM.PixelRegion); %read subsetted DEM

    % load, subset dH
    DH.info=geotiffinfo(DH.path);
    DH.R = DH.info.RefMatrix;
    [DH.PixelRegion,DH.LatG,DH.LonG] = subset_geo(DH.info,BBoxLL);
    DH.data=imread(DH.path,'PixelRegion',DH.PixelRegion);
    
    if dhfilter==1 %if noisy or gappy, apply some preprocessing to the dH
        DH.data2=DH.data;
        DH.errthresh1=3.*nanstd(DH.data(:));
    %     DH.MM=movmean(DH.data,[9 9],'omitnan'); %calc local mean
    %     DH.MST=movmean(DH.data,[9 9],'omitnan'); %calc local std
    %     DH.MSG = (DH.data-DH.MM)./DH.MST; %stdev from local mean
    %     DH.data2(abs(DH.MSG)>2)=NaN; %filter sigma>2
    %     DH.dH3=inpaint_nans(DH.data2);
        DH.data2(abs(DH.data2)>DH.errthresh1)=NaN;
        DH.errthresh2=3.*nanstd(DH.data2(:));
        DH.data2(abs(DH.data2)>DH.errthresh2)=NaN;
        DH.dH3=imgaussfilt(DH.data2); %smooths dH slightly, expands NaN around bad DH data
    else
        DH.dH3=DH.data;
    end
    
    %% REPROJECT ALL (to THX coordinate system)

    buffdist = 1000; %expands the domain for subsetting around the glacier

    %new coordinates based on the THX data, but at DX interval
    N.x3 = DX.*[(floor((THX.xm(1)-buffdist)/DX)):((ceil(THX.xm(end)+buffdist)/DX))];
    N.y3 = DX.*[(ceil((THX.ym(1)+buffdist)/DX)):-1:(floor((THX.ym(end)-buffdist)/DX))];
    [N.x3g,N.y3g] = meshgrid(N.x3,N.y3);

    N.Rout = [0,-DX;DX,0;N.x3(1)-DX,N.y3(1)+DX]; %matlab georeferencing is off

    %RESAMPLE THICKNESS
%     [THX.xN,THX.yN] = projfwd(THX.info,THX.LatG,THX.LonG); %compute projected coordinates
%     N.THX = griddata(THX.xN(:),THX.yN(:),double(THX.data(:)),N.x3g(:),N.y3g(:),'cubic');
    N.THX = griddata(THX.xmG(:),THX.ymG(:),double(THX.data(:)),N.x3g(:),N.y3g(:),'cubic');
    N.THX = reshape(N.THX,size(N.x3g));

    %resample velocity data
    [V.xN,V.yN] = projfwd(THX.info,V.LatG(:),V.LonG(:)); %compute projected coordinates
    V.xN(iERR)=[];V.yN(iERR)=[];
    if isempty(V.info.PCS) %if velocity data has a source projection
        V.Uraw2=V.Uraw(:);V.Uraw2(iERR)=[]; %OLD
        V.Vraw2=V.Vraw(:);V.Vraw2(iERR)=[]; %OLD
    else %geogrpahic - assume values are relative to N and projected!!
        [V.x2N,V.y2N] = projfwd(THX.info,V.Lat2(:),V.Lon2(:)); %compute projected coordinates for vector end-points
        V.x2N(iERR)=[];V.y2N(iERR)=[];
        V.Uraw2=V.x2N-V.xN; %velocity vector in new coord system
        V.Vraw2=V.y2N-V.yN; %velocity vector in new coord system
    end
    N.U = griddata(V.xN(:),V.yN(:),double(V.Uraw2(:)),N.x3g(:),N.y3g(:),'cubic');
    N.U = reshape(N.U,size(N.x3g));
    N.V = griddata(V.xN(:),V.yN(:),double(V.Vraw2(:)),N.x3g(:),N.y3g(:),'cubic');
    N.V = reshape(N.V,size(N.x3g));
    N.S = sqrt(N.V.^2+N.U.^2); %speed as velocity vector magnitude

    %resample DEM
    [DEM.xN,DEM.yN] = projfwd(THX.info,DEM.LatG,DEM.LonG); %compute projected coordinates
    N.DEM = griddata(DEM.xN(:),DEM.yN(:),double(DEM.data(:)),N.x3g(:),N.y3g(:),'cubic');
    N.DEM = reshape(N.DEM,size(N.x3g));

    %resample dH
    [DH.xN,DH.yN] = projfwd(THX.info,DH.LatG,DH.LonG); %compute projected coordinates
    N.DH = griddata(DH.xN(:),DH.yN(:),double(DH.dH3(:)),N.x3g(:),N.y3g(:),'cubic');
    N.DH = reshape(N.DH,size(N.x3g));

    %set mask
    N.MASK = N.THX>0;

    %gap-fill dH based on elevation
    N.DH0=N.DH; %unfilled values
    N.DH(N.MASK)= index_nanfill(N.DH(N.MASK),N.DEM(N.MASK));
    
    % MASK
    N.DH((N.MASK==0))=0;
    N.U((N.MASK==0))=0;
    N.V((N.MASK==0))=0;
    N.THX((N.MASK==0))=0;

    %% FLUX AND EMERGENCE CALCS - A LOGICAL PLACE TO LOOP ITERATIONS?
    % convert velocity to column-averaged with a fixed multiplier - 0.9 for now but can be optimized
    umult=0.9; %range is [0.8,1], more realistically [0.81,0.99]
    N.Smean=umult.*N.S;
    N.Umean=umult.*N.U;
    N.Vmean=umult.*N.V;
    
    %% calculate flux divergence per pixel
    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));
        
    %initialize
    N.FDIV = zeros(size(N.THX));
    N.FDIVx = zeros(size(N.THX));
    N.FDIVy = zeros(size(N.THX));

    %pixel-based flux magnitude
    N.FLUX = N.Smean.*N.THX;N.FLUX(N.FLUX<=0)=NaN;
    
    %first order centered-difference
    N.FDIVx(:,2:end-1) = (N.Umean(:,3:end).*N.THX(:,3:end)-N.Umean(:,1:end-2).*N.THX(:,1:end-2))/2./dx;
    N.FDIVy(2:end-1,:) = -(N.Vmean(3:end,:).*N.THX(3:end,:)-N.Vmean(1:end-2,:).*N.THX(1:end-2,:))/2./dx; %negative because V is opposite direction to raster pixels
    N.FDIV = (N.FDIVx+N.FDIVy);%total and normalize to area, m/yr

    %trim to mask
    N.FDIV(N.MASK==0)=NaN;

    cd(outdir)

    %% zone segmentation, calc fluxes
    N.zones = uint16(segment_Gmask_slope2(N.DEM,N.MASK,DX,segdist)); %segments glacier mask and DEM into approximately uniformly-spaced elevation bands

    %aggregate variables over zones
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
    ind4 = N.DEM>median(N.DEM(N.MASK));
    
    N.Hdensity(ind1&~ind2)=0.9; %thinning and emergence = melt
    N.Hdensity(~ind1&ind2)=0.6; %thickening and submergence = acc
    N.Hdensity(ind1&ind2&ind3)=0.6; %emergence and thickening, more thickening - acc
    N.Hdensity(ind1&ind2&~ind3)=0.85; %emergence and thickening, more thickening - mixed
    N.Hdensity(~ind1&~ind2&ind3)=0.9; %submergence and thinning, more thinning - melt
    N.Hdensity(~ind1&~ind2&~ind3)=0.85; %submergence and thinning, less thinning - mixed
    N.Hdensity((ind4==0)&N.MASK)=0.9; %
   
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
    
%%
if plotouts==1
    cmap = cbrewer('div','RdBu',21);
    
    cmap2 = [flipud(cbrewer('seq','Reds',11));cbrewer('seq','Blues',5)];
% cmap = [0,0,0;cbrewer('div','RdBu',21);0,0,0];

    Nout = 'grid-emergence'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(-1.*N.FDIV)
    colorbar;colormap(flipud(cmap))
    title(Ntitle)
    caxis([-10,10])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    %
    Nout = 'EL-zones'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.zones)
    colorbar;colormap([0,0,0;lines(256)])
    title(Ntitle)
%     caxis([-5,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    
    %
    Nout = 'EL-zone-avg-emergence'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(-1.*N.z2fdiv)
    colorbar;colormap(flipud(cmap))
    title(Ntitle)
    caxis([-5,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
   
    Nout = 'grid-dh'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.DH)
    colorbar;colormap(cmap2)
    title(Ntitle)
    caxis([-10,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])

    Nout = 'grid-SMB'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.SMB)
    colorbar;colormap(cmap2)
    title(Ntitle)
    caxis([-10,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    %
    Nout = 'grid-Hdensity'
    Ntitle = [Nout ' (1000 km m^{-3}), ' outtitle1];
    figure
    imagesc((N.Hdensity))
    colorbar;
    title(Ntitle)
    caxis([.5,1])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    %
    Nout = 'thx'
    Ntitle = [Nout ' (m), ' outtitle1];
    figure
    imagesc(N.THX)
    colorbar;%colormap(cmap)
    title(Ntitle)
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])

    %
    Nout = 'EL-zonal-SMB'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.SMBz2)
    colorbar;colormap(cmap2)
    title(Ntitle)
    caxis([-10,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])

    %
    Nout = 'surface-speed'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.Smean)
    colorbar;
    title(Ntitle)
    caxis([0,100])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])

end
    %% plot SMB curves and fluxes with elevation

    [ELA,FLout,cFLu,tFL,ELs,ELfluxes]=through_fluxes(N.MASK,N.DEM,N.Umean,N.Vmean,N.THX,DX);%,dy,sig_H,UE,VE,ERRORs)
    
    %
    Nout = 'SMB-EL'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    iN = find(N.MASK);
    figure;
    scatter(N.DEM(iN),N.SMB(iN),4,'kd');hold on
    scatter(N.DEM(iN),N.SMBz2(iN),4,'b*');
    ylabel('SMB m/a')
    ylim([-20,10])    
    title(Ntitle)
    legend('gridded SMB','zonal SMB','elevation fluxes')
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
%%
    dZ=25;
    ELS=dZ*[floor(min(N.DEM(iN)/dZ)):ceil(max(N.DEM(iN)/dZ))];
    elSMBg = zeros(1,length(ELS)-1);
    elSMBz = zeros(1,length(ELS)-1);
    elSMBe = zeros(1,length(ELS)-1);
    hyps = zeros(1,length(ELS)-1);
    for iel=1:length(ELS)-1
        cIN=N.MASK&(N.DEM>=ELS(iel))&(N.DEM<ELS(iel+1));
        elSMBg(iel)=nanmean(N.SMB(cIN));
        elSMBz(iel)=nanmean(N.SMBz2(cIN));
        elSMBe(iel)=nanstd(N.SMB(cIN));
    end
    ELS2=unique(dZ.*floor(N.DEM(iN)./dZ));
    
    Nout = 'SMB-EL-fluxes'
    Ntitle = [Nout ', ' outtitle1];
    figure;
    yyaxis left
%     plot(ELS(1:end-1)+dZ/2,elSMBg,'b');hold on
%     plot(ELS(1:end-1)+dZ/2,elSMBz,'k');hold on
    boxplot(N.SMB(iN),dZ.*floor(N.DEM(iN)./dZ),'Positions',ELS2);hold on
    ylabel('SMB m/a')
    ylim([-20,10])
    yyaxis right
    plot(ELs,ELfluxes,'linewidth',2)
    ylim([0,3*max(ELfluxes)])
    ylabel('Flux m^3/a')
    title(Ntitle)
%     legend('gridded','zonal','flux','location','southeast')
%     legend('flux','location','southeast')
    set(gca,'XTick',1000.*floor(min(N.DEM(iN)./1000)):500:1000.*ceil(max(N.DEM(iN)./1000)))
    set(gca,'XTickLabel',get(gca,'XTick'))
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    
    %% export 
if exports==1

    % DEM
    geotiffwrite([Glacier '_DEM.tif'],uint16(N.DEM),N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % THX
    geotiffwrite([Glacier '_THX.tif'],N.THX,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % VEL
    geotiffwrite([Glacier '_Smean.tif'],N.Smean,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % dH
    geotiffwrite([Glacier '_dH.tif'],N.DH,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % FDIV
    geotiffwrite([Glacier '_FDIV.tif'],N.FDIV,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % SMB
    geotiffwrite([Glacier '_SMB.tif'],N.SMB,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % Zones
    geotiffwrite([Glacier '_zones.tif'],uint16(N.zones),N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % zFDIV
    geotiffwrite([Glacier '_zFDIV.tif'],N.z2fdiv,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % zSMB
    geotiffwrite([Glacier '_zSMB.tif'],N.SMBz2,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

end

%%    

disp('finished')
   
