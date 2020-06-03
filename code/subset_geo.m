function [PixelRegion,LatG,LonG,xmG1,ymG1] = subset_geo(info,BBoxLL)
%subset_geo - Function to determine pixel and map coordinates for a raster
%within a bounding box of geographic coordinates
%
% Syntax:  [PixelRegion,LatG,LonG,xmG1,ymG1] = subset_geo(info,BBoxLL)
%
% Inputs:
%    info - geotiffinfo structure for raster to be subset
%    BBoxLL - 2x2 matrix [Latmin,Latmax;Lonmin,Lonmax] of an ROI
%
% Outputs:
%    PixelRegion - pixel indices for raster for rectangle containing bbox [y1,yend;x1,xend]
%    LatG - Latitude values for PixelRegion
%    LonG - Longitude values for PixelRegion
%    xmG1 - x map coordinates for PixelRegion
%    ymG1 - y map coordinates for PixelRegion
%
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-June-2020

%% ------------- BEGIN CODE --------------
%pixel map coordinates for raster to be subset
xp = [1:info.Width];
yp = [1:info.Height]; %pixel coordinates are cell centers
[xm,~] = pix2map(info.RefMatrix,ones(size(xp)),xp);
[~,ym] = pix2map(info.RefMatrix,yp,ones(size(yp)));
[xmG,ymG] = meshgrid(xm,ym);

%transform to geographic coordinates 
if isempty(info.GeoTIFFCodes.PCS) %already a GCS
    LatG1=ymG;
    LonG1=xmG;
else
    [LatG1,LonG1] = projinv(info,xmG(:),ymG(:)); %inverse of projection to WGS84
end

%column returned, reshape to grid
LatG1=reshape(LatG1,size(xmG));
LonG1=reshape(LonG1,size(xmG));

%find area within BBox
LatReg = max((LatG1>BBoxLL(3))&(LatG1<BBoxLL(4)),[],2);
LonReg = max((LonG1>BBoxLL(1))&(LonG1<BBoxLL(2)),[],1);

%find indices
y1=find(LatReg,1,'first');
ye=find(LatReg,1,'last');
x1=find(LonReg,1,'first');
xe=find(LonReg,1,'last');
PixelRegion = {[y1,ye],[x1,xe]};
m=ye-y1+1;
n=xe-x1+1;

%grid of Lat and Lon for subset area
LatG=reshape(LatG1(LatReg&LonReg),m,n);
LonG=reshape(LonG1(LatReg&LonReg),m,n);

%grid of map coordinates for subset area
xmG1=xmG(y1:ye,x1:xe);
ymG1=ymG(y1:ye,x1:xe);
