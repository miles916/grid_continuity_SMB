clear 
close all

homedir = '\\wsl.ch\fe\gebhyd\8_Him\Personal_folders\Evan\flux_thickness_iteration'; %EACH USER TO UPDATE
addpath(genpath([homedir '\code']))

plotouts=1; %output plots or not
exports=1; %save geotiffs or not
DX = 100; %resolution to run calculations at 
Glacier = 'Rongbuk'
% Glacier = 'Matanuska'
datatitle = ['test_' Glacier];


outdir = [homedir '\results\' datatitle '_' date]
mkdir(outdir)

%% configure inputs
switch Glacier
    case 'Rongbuk'
        V.pathx = fullfile(homedir,'data',Glacier,'HMA_G0120_vx.tif');
%         V.pathxerr = fullpath(homedir,'data',Glacier,'HMA_G0120_vx_err.tif');
        V.pathy = fullfile(homedir,'data',Glacier,'HMA_G0120_vy.tif');
%         V.pathyerr = fullpath(homedir,'data',Glacier,'HMA_G0120_vy_err.tif');
        DH.path = fullfile(homedir,'data',Glacier,'15.09991_dH.tif');
%         THX.path = fullfile(homedir,'data',Glacier,'RGI60-15.09991_thickness_composite.tif');
        THX.path = fullfile(homedir,'data',Glacier,'thickness_RGI60-15.09991_HF2012.tif');
        DEM.path = fullfile(homedir,'data',Glacier,'15.09991_AW3D.tif');
        segdist = 300;
        V.Vmult=1; % velocity axes in correct direction
    case 'Matanuska'
        V.pathx = fullfile(homedir,'data',Glacier,'Vx_Matanuska_ITSLIVE.tif');
        V.pathy = fullfile(homedir,'data',Glacier,'Vy_Matanuska_ITSLIVE.tif');
        DH.path = fullfile(homedir,'data',Glacier,'dhdt_Matanuska_2007_2016.tif');
        THX.path = fullfile(homedir,'data',Glacier,'RGI60-01.10557_thickness_composite.tif');
        THX.path = fullfile(homedir,'data',Glacier,'thickness_RGI60-01.10557_HF2012.tif');
        DEM.path = fullfile(homedir,'data',Glacier,'GDEM_Matanuska.tif');
        segdist=2000;
        V.Vmult=1; % velocity y-axis reversed to grid direction
end


    outtitle1 = [num2str(DX) 'mgrid'];
    cd(outdir)
    
    %% Load all input data, resample to common grid, etc
    % load THX, determine LL bounding box
    THX.data=geotiffread(THX.path);
    THX.info=geotiffinfo(THX.path);
    THX.R = THX.info.RefMatrix;
    THX.xp=1:THX.info.Width;
    THX.yp=1:THX.info.Height;
    [THX.xm,~] = pix2map(THX.R,ones(size(THX.xp)),THX.xp);
    [~,THX.ym] = pix2map(THX.R,THX.yp,ones(size(THX.yp)));
    [THX.xmG,THX.ymG] = meshgrid(THX.xm,THX.ym);
    [THX.LatG,THX.LonG] = projinv(THX.info,THX.xmG(:),THX.ymG(:));
    
    % set initial mask from thickness data
    MASK0 = THX.data>0;
    
    BBoxUTM = THX.info.BoundingBox;
    BBoxLL = [min(THX.LonG(MASK0)) min(THX.LatG(MASK0));max(THX.LonG(MASK0)) max(THX.LatG(MASK0))];
    
    % load, subset velocity
    V.info = geotiffinfo(V.pathx); %easiest way to get correct projection details
    V.R = V.info.RefMatrix;
    V.xp=1:V.info.Width;
    V.yp=1:V.info.Height;
    [V.xm,~] = pix2map(V.R,ones(size(V.xp)),V.xp);
    [~,V.ym] = pix2map(V.R,V.yp,ones(size(V.yp)));
    [V.PixelRegion,V.LatG,V.LonG,V.xmG,V.ymG] = subset_geo(V.info,BBoxLL);
    V.Uraw=imread(V.pathx,'PixelRegion',V.PixelRegion);
    V.Vraw=imread(V.pathy,'PixelRegion',V.PixelRegion).*V.Vmult;
    
    %calculate velocity end-points
    V.xm2=V.xmG+double(V.Uraw); %velocity vector end-point
    V.ym2=V.ymG+double(V.Vraw); %velocity vector end-point
    [V.Lat2,V.Lon2]=projinv(V.info,V.xm2(:),V.ym2(:)); %unprojected velocity vector end-point
    
    %identify likely errors
    iERR=(abs(V.Uraw)>400)|(abs(V.Vraw)>400);
    
%     histogram(V.Vraw(iERR==0),100)
    
    % load, subset DEM
    try
        DEM.info=geotiffinfo(DEM.path);
    catch
        DEM.info=imfinfo(DEM.path);
        DEM.info=DEM.info(1);
        DEM.info.RefMatrix = [0,-DEM.info.ModelPixelScaleTag(1);DEM.info.ModelPixelScaleTag(2),0;DEM.info.ModelTiepointTag(4),DEM.info.ModelTiepointTag(5)];
    end
    DEM.R = DEM.info.RefMatrix;
    [DEM.PixelRegion,DEM.LatG,DEM.LonG] = subset_geo(DEM.info,BBoxLL);
    DEM.data=imread(DEM.path,'PixelRegion',DEM.PixelRegion);

    % load, subset dH
    DH.info=geotiffinfo(DH.path);
    DH.R = DH.info.RefMatrix;
    [DH.PixelRegion,DH.LatG,DH.LonG] = subset_geo(DH.info,BBoxLL);
    DH.data=imread(DH.path,'PixelRegion',DH.PixelRegion);
    
    DH.data2=DH.data;
    DH.MM=movmean(DH.data,[9 9],'omitnan'); %calc local mean
    DH.MST=movmean(DH.data,[9 9],'omitnan'); %calc local std
    DH.MSG = (DH.data-DH.MM)./DH.MST; %stdev from local mean
    DH.data2(abs(DH.MSG)>2)=NaN; %filter sigma>2
    DH.dH3=inpaint_nans(DH.data2);

    %% REPROJECT ALL (to THX coordinate system)

    buffdist = 1000;
    %new coordinates based on the THX data, but at DX interval
    N.x3 = DX.*[(floor((THX.xm(1)-buffdist)/DX)):((ceil(THX.xm(end)+buffdist)/DX))];
    N.y3 = DX.*[(ceil((THX.ym(1)+buffdist)/DX)):-1:(floor((THX.ym(end)-buffdist)/DX))];
    [N.x3g,N.y3g] = meshgrid(N.x3,N.y3);

    N.Rout = [0,-DX;DX,0;N.x3(1),N.y3(1)];

    %RESAMPLE THICKNESS
%     [THX.xN,THX.yN] = projfwd(THX.info,THX.LatG,THX.LonG); %compute projected coordinates
%     N.THX = griddata(THX.xN(:),THX.yN(:),double(THX.data(:)),N.x3g(:),N.y3g(:),'cubic');
    N.THX = griddata(THX.xmG(:),THX.ymG(:),double(THX.data(:)),N.x3g(:),N.y3g(:),'cubic');
    N.THX = reshape(N.THX,size(N.x3g));

    %resample velocity data
    [V.xN,V.yN] = projfwd(THX.info,V.LatG(:),V.LonG(:)); %compute projected coordinates
    [V.x2N,V.y2N] = projfwd(THX.info,V.Lat2(:),V.Lon2(:)); %compute projected coordinates for vector end-points
    V.xN(iERR)=[];V.yN(iERR)=[];
    V.x2N(iERR)=[];V.y2N(iERR)=[];
%     V.Uraw2=V.Uraw(:);V.Uraw2(iERR)=[]; %OLD
%     V.Vraw2=V.Vraw(:);V.Vraw2(iERR)=[]; %OLD
    V.Uraw2=V.x2N-V.xN; %velocity vector in new coord system
    V.Vraw2=V.y2N-V.yN; %velocity vector in new coord system
    N.U = griddata(V.xN(:),V.yN(:),double(V.Uraw2(:)),N.x3g(:),N.y3g(:),'cubic');
    N.U = reshape(N.U,size(N.x3g));
    N.V = griddata(V.xN(:),V.yN(:),double(V.Vraw2(:)),N.x3g(:),N.y3g(:),'cubic');
    N.V = reshape(N.V,size(N.x3g));
    N.S = sqrt(N.V.^2+N.U.^2);

    %resample DEM
    [DEM.xN,DEM.yN] = projfwd(THX.info,DEM.LatG,DEM.LonG); %compute projected coordinates
    N.DEM = griddata(DEM.xN(:),DEM.yN(:),double(DEM.data(:)),N.x3g(:),N.y3g(:),'cubic');
    N.DEM = reshape(N.DEM,size(N.x3g));

    %resample dH
    [DH.xN,DH.yN] = projfwd(THX.info,DH.LatG,DH.LonG); %compute projected coordinates
    N.DH = griddata(DH.xN(:),DH.yN(:),double(DH.dH3(:)),N.x3g(:),N.y3g(:),'cubic');
    N.DH = reshape(N.DH,size(N.x3g));

    N.MASK = N.THX>0;

    % MASK

    N.DH((N.MASK==0))=0;
    N.U((N.MASK==0))=0;
    N.V((N.MASK==0))=0;
    N.THX((N.MASK==0))=0;

    %% FLUX AND EMERGENCE CALCS - A LOGICAL PLACE TO LOOP ITERATIONS?
    % convert velocity to column-averaged
   
    %assume fixed multiplier
    umult=0.9; %range is [0.8,1], more realistically [0.81,0.99], (uncertainty is ~0.1)
    N.Smean=umult.*N.S;
    N.Umean=umult.*N.U;
    N.Vmean=umult.*N.V;
%     sig_umult = 0.1;
    
    %% calculate flux divergence per pixel
    dx = mode(diff(N.x3));
    dy = mode(diff(N.y3));
        
    N.FDIV = zeros(size(N.THX));
    N.FDIVx = zeros(size(N.THX));
    N.FDIVy = zeros(size(N.THX));

    N.FLUX = N.Smean.*N.THX;N.FLUX(N.FLUX<=0)=NaN;
    
    %first order centered-difference
    N.FDIVx(:,2:end-1) = (N.Umean(:,1:end-2).*N.THX(:,1:end-2)-N.Umean(:,3:end).*N.THX(:,3:end)).*(dy)/2;
    N.FDIVy(2:end-1,:) = (N.Vmean(1:end-2,:).*N.THX(1:end-2,:)-N.Vmean(3:end,:).*N.THX(3:end,:)).*(dx)/2;
    N.FDIV = (N.FDIVx+N.FDIVy)/(dy*dx);%total and normalize to area, m/yr

    N.FDIV(N.MASK==0)=NaN;

    cd(outdir)

    %% zone segmentation, calc fluxes
    N.zones = uint16(segment_Gmask_slope2(N.DEM,N.MASK,DX,segdist));

    %aggregate variables over zones
    N.z2fdiv = zonal_aggregate(N.zones,N.FDIV); %simply aggregates values in the zone - FAST simple mean; not quite robust but produces the same value as a perimeter integration (good...)
    N.z2DH = zonal_aggregate(N.zones,N.DH); %simply aggregates values in the zone - FAST

     %% SMB and uncertainty
    
    N.SMB = N.DH-N.FDIV;
    N.SMB(N.MASK==0)=NaN;

    N.SMBz2 = N.z2DH-N.z2fdiv;


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
% cmap = [0,0,0;cbrewer('div','RdBu',21);0,0,0];

    Nout = 'grid-emergence'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.FDIV)
    colorbar;colormap(flipud(cmap))
    title(Ntitle)
    caxis([-5,5])
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
    imagesc(N.z2fdiv)
    colorbar;colormap(flipud(cmap))
    title(Ntitle)
    caxis([-5,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
   
    Nout = 'grid-dh'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.DH)
    colorbar;colormap(cmap)
    title(Ntitle)
    caxis([-5,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])

    Nout = 'grid-SMB'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.SMB)
    colorbar;colormap(cmap)
    title(Ntitle)
    caxis([-5,5])
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
    colorbar;colormap(cmap)
    title(Ntitle)
    caxis([-5,5])
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
    legend('hybrid zonal SMB','zonal SMB','elevation fluxes')
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
%%
    dZ=10;
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
    
    Nout = 'SMB-EL-fluxes'
    Ntitle = [Nout ', ' outtitle1];
    figure;
    yyaxis left
    plot(ELS(1:end-1)+dZ/2,elSMBg,'b');hold on
    plot(ELS(1:end-1)+dZ/2,elSMBz,'k');hold on
    ylabel('SMB m/a')
    yyaxis right
    plot(ELs,ELfluxes,'linewidth',2)
    ylim([0,3*max(ELfluxes)])
    ylabel('Flux m^3/a')
    title(Ntitle)
    legend('gridded','hybrid','flux','location','southeast')
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    
    %% export 
if exports==1

    % DEM
    geotiffwrite([Glacier '_DEM.tif'],uint16(N.DEM),N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % THX
    geotiffwrite([Glacier '_THX.tif'],N.THX,N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % VEL
    geotiffwrite([Glacier '_Smean.tif'],N.Smean,N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % dH
    geotiffwrite([Glacier '_dH.tif'],N.DH,N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % FDIV
    geotiffwrite([Glacier '_FDIV.tif'],N.FDIV,N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)])

    % SMB
    geotiffwrite([Glacier '_SMB.tif'],N.SMB,N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)])

    % Zones
    geotiffwrite([Glacier '_zones.tif'],uint16(N.zones),N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % zFDIV
    geotiffwrite([Glacier '_zFDIV.tif'],N.z2fdiv,N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % zSMB
    geotiffwrite([Glacier '_zSMB.tif'],N.SMBz2,N.Rout,'CoordRefSysCode',['EPSG:' num2str(THX.info.GeoTIFFCodes.PCS)],'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

end

%%    

disp('finished')
   
