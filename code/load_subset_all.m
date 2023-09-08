function [THX,DEM,V,DH]=load_subset_all(P,THX,DEM,V,DH)
    % load THX, determine LL bounding box
    THX.data=geotiffread(THX.path); %read data
    THX.info=geotiffinfo(THX.path); %read metadata
    THX.R = THX.info.RefMatrix; %georeferencing matrix
    
        %setup grid of coordinates for thickness data
%     THX.xp=1:THX.info.Width;
%     THX.yp=1:THX.info.Height;
%     [THX.xm,~] = pix2map(THX.R,ones(size(THX.xp)),THX.xp);
%     [~,THX.ym] = pix2map(THX.R,THX.yp,ones(size(THX.yp)));
    [THX.xm,THX.ym] = pixcenters(THX.R,[THX.info.Height,THX.info.Width]);
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
    if isfield(V,'pathxe')
        V.UrawE=V.mult.*imread(V.pathxe,'PixelRegion',V.PixelRegion);
        V.VrawE=V.mult.*imread(V.pathye,'PixelRegion',V.PixelRegion);
    end
    
    V.Vreproj=P.Vreproj;
    if P.Vreproj == 1
            %calculate velocity end-points (reproject full velocity vector)
        V.xm2=V.xmG+double(V.Uraw); %velocity vector end-point
        V.ym2=V.ymG+double(V.Vraw); %velocity vector end-point
        if isempty(V.info.PCS)==0 %if velocity product is projected, reproject the vectors. if not, we assume it to be oriented north and with units of m/a already
            [V.Lat2,V.Lon2]=projinv(V.info,V.xm2(:),V.ym2(:)); %unprojected velocity vector end-point
        end
    end
    
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
    DH.data=imread(DH.path,'PixelRegion',DH.PixelRegion).*DH.mult;
