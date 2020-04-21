function [PixelRegion,LatG,LonG,xmG1,ymG1] = subset_geo(info,BBoxLL)

xp = [1:info.Width];
yp = [1:info.Height]; %pixel coordinates are cell centers
[xm,~] = pix2map(info.RefMatrix,ones(size(xp)),xp);
[~,ym] = pix2map(info.RefMatrix,yp,ones(size(yp)));
[xmG,ymG] = meshgrid(xm,ym);
if isempty(info.GeoTIFFCodes.PCS)
    LatG1=ymG;
    LonG1=xmG;
else
    [LatG1,LonG1] = projinv(info,xmG(:),ymG(:));
end
    
   
LatG1=reshape(LatG1,size(xmG));
LonG1=reshape(LonG1,size(xmG));

LatReg = max((LatG1>BBoxLL(3))&(LatG1<BBoxLL(4)),[],2);
LonReg = max((LonG1>BBoxLL(1))&(LonG1<BBoxLL(2)),[],1);

y1=find(LatReg,1,'first');
ye=find(LatReg,1,'last');
x1=find(LonReg,1,'first');
xe=find(LonReg,1,'last');
PixelRegion = {[y1,ye],[x1,xe]};
m=ye-y1+1;
n=xe-x1+1;

LatG=reshape(LatG1(LatReg&LonReg),m,n);
LonG=reshape(LonG1(LatReg&LonReg),m,n);

xmG1=xmG(y1:ye,x1:xe);
ymG1=ymG(y1:ye,x1:xe);
