%master_ContinuitySMB_ex.m - Master script to estimate surface mass balance
%distribution for a glacier from inputs of ice thickness, thinning, and
%velocity, based on the continuity equation (see, e.g. Bisset et al, 2020: https://doi.org/10.3390/rs12101563)
%
% New glaciers require this template to be coded with paths for the full set of
% inputs, and any additional needed preprocessing. 
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

%glacier ID
Glacier = 'Rongbuk'

%title of output directory
datatitle = ['test_' Glacier];

%inputs        
V.pathx = fullfile(homedir,'data',Glacier,'HMA_G0120_vx.tif');
V.pathy = fullfile(homedir,'data',Glacier,'HMA_G0120_vy.tif');
DH.path = fullfile(homedir,'data',Glacier,'15.09991_dH.tif');
THX.path = fullfile(homedir,'data',Glacier,'thickness_RGI60-15.09991_HF2012.tif');
DEM.path = fullfile(homedir,'data',Glacier,'15.09991_AW3D.tif');

%settings
plotouts=1; %output plots or not
exports=1; %save geotiffs or not
DX = 200; %resolution to run calculations at 
segdist=300; %effective linear distance between flowbands
V.mult=1; %scale to convert input units to m/a
V.filter=0; %swtich to smooth velocity data or not
dhfilter=0; %switch to peform 2x 3-sigma outlier removal (from overall dataset - only if erroneous pixels are common)
umult=2; %switch for column-average velocity [0.8-1] is physical range. '0' estimates based on THX dist. '2' estimates for each pixel. 


%% initialization
    addpath(genpath([homedir '\code'])) %add path to related scripts

    outdir = [homedir '\results\' datatitle '_' date]
    mkdir(outdir)

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
    if V.filter==1 %gaussian low-pass filter removing extreme variations
        V.Uraw=imgaussfilt(V.Uraw,5);
        V.Vraw=imgaussfilt(V.Vraw,5);
    end
    
        %calculate velocity end-points (reproject full velocity vector)
    V.xm2=V.xmG+double(V.Uraw); %velocity vector end-point
    V.ym2=V.ymG+double(V.Vraw); %velocity vector end-point
    if isempty(V.info.PCS)==0 %if velocity product is projected, reproject the vectors. if not, we assume it to be oriented north and with units of m/a already
        [V.Lat2,V.Lon2]=projinv(V.info,V.xm2(:),V.ym2(:)); %unprojected velocity vector end-point
    end
    
    %identify likely errors
    V.iERR=(abs(V.Uraw)>400)|(abs(V.Vraw)>400);
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
    
    %% Cleaning raw data if needed
    
        % filter thickness data
    THX.dataR=THX.data;
    THX.data=imgaussfilt(THX.dataR,4); %gaussian low-pass filter. Important for thickness maps derived from field data
    THX.data(MASK0==0)=0;
    
    if dhfilter==1 %if noisy or gappy, apply some preprocessing to the dH
        DH.data2=DH.data;
        DH.errthresh1=3.*nanstd(DH.data(:));
        DH.data2(abs(DH.data2)>DH.errthresh1)=NaN;
        DH.errthresh2=3.*nanstd(DH.data2(:));
        DH.data2(abs(DH.data2)>DH.errthresh2)=NaN;
        DH.dH3=imgaussfilt(DH.data2); %smooths dH slightly, expands NaN around bad DH data
    else
        DH.dH3=DH.data;
    end
    
    %% REPROJECT AND RESAMPLE ALL INPUTS (to THX coordinate system)
    N = resample_inputs(DX,THX,V,DEM,DH);

    %set mask
    N.MASK = N.THX>0;

    % MASK
    N.DH((N.MASK==0))=0;
    N.U((N.MASK==0))=0;
    N.V((N.MASK==0))=0;
    N.THX((N.MASK==0))=0;
    
    %% clean reprojected inputs if needed
    
    %gap-fill dH based on elevation
    N.DH0=N.DH; %unfilled values
    N.DH(N.MASK)= index_nanfill(N.DH(N.MASK),N.DEM(N.MASK));
    N.DH((N.MASK==0))=0;

    %% zone segmentation
    N.zones = uint16(segment_Gmask_slope2(N.DEM,N.MASK,N.DX,segdist)); %segments glacier mask and DEM into approximately uniformly-spaced elevation bands

    %% determine column-average velocity correction
    % convert velocity to column-averaged with a fixed multiplier - 0.9 for now but can be optimized
    if umult==0 %estimates based on range of ice thicknesses
        [umult,~,~]=f_probability2(N.THX(N.MASK)); %optimizes based on range of ice thicknesses
%         umult=0.9; %range is [0.8,1], more realistically [0.81,0.99]
    elseif umult==2 %estimates based ice thickness of each individual pixel
        [umult,~]=f_probability1(N.THX,N.MASK); %optimizes based on each ice thickness
    elseif umult<0.8|umult>1 
        error('Unphyiscal multiplier for column-averaged velocity')
    end

    %% SMB calculations
    cd(outdir)
    N=FluxCalcsSimple(N,umult); %calculates fluxes 
    
%% DETERMINE FLUXES THROUGH EACH ELEVATION BAND
    [N.ELA,N.FLout,N.cFLu,N.tFL,N.ELs,N.ELfluxes]=through_fluxes(N.MASK,N.DEM,N.Umean,N.Vmean,N.THX,N.DX);%,dy,sig_H,UE,VE,ERRORs)

if plotouts==1
    SMBplots(N,outtitle1,Glacier);
end
    
    %% export geotifs
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
   
