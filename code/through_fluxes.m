function [ELA,FLout,cFLu,tFL,ELs,Afluxes]=through_fluxes(mask,DEM,Umean,Vmean,thx,dx)
% through_fluxes - Function to calculate volumetric fluxes of ice with 
% respect to elevation. Corrects x- and y- components of velocity into 
% down-and cross-glacier components based on smoothed DEM surface gradient,
% then integrates down-glacier component across 25m elevation contours
%
% Syntax:  [ELA,FLout,cFLu,tFL,ELs,Afluxes]=through_fluxes(mask,DEM,Umean,Vmean,thx,dx)
%
% Inputs:
%    mask - binary raster mask of glacier extent
%    DEM - double-precision raster of elevation with same extent as mask
%    Umean - double-precision raster of x-component of column-averaged velocity with same extent as mask
%    Vmean - double-precision raster of y-component of column-averaged velocity with same extent as mask
%    thx - double-precision raster of ice thickness with same extent as mask
%    dx - pixel size in porjected units corresponding to veloctiy and DEM units (m)
%
% Outputs:
%    ELA - elevation of peak down-glacier flux
%    FLout - peak down-glacier ice flux (m3/a) 
%    cFLu - raster of down-glacier oriented volumetric fluxes
%    tFL - raster of volumetric flux at each pixel
%    ELs - array of elevations at 25m intervals
%    Afluxes - array of volumetric ice fluxes through ELs
%
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 01-July-2020

DEM1=DEM;

DEM1(mask==0)=NaN;

DEM2 = inpaint_nans(DEM1,0); %repeats edges values

DEM3 = imgaussfilt(DEM2,2);
% DEM3 = imgaussfilt(DEM3,2);

ELs = 25.*(ceil(min(DEM3(DEM3(:).*mask(:)>0))/25):floor(max(DEM3(:).*mask(:))/25));

[Fx,Fy] = gradient(DEM3);
SLO = sqrt(Fx.^2+Fy.^2)/dx;
SLO(mask==0)=NaN;
Fx(mask==0)=NaN;
Fy(mask==0)=NaN;
DEM3(mask==0)=NaN;
fx = Fx./sqrt(Fx.^2+Fy.^2);
fy = Fy./sqrt(Fx.^2+Fy.^2);

%% calculate flux corrected for slope direction 
%     FLx=Umean.*thx;
%     FLy=Vmean.*thx;
    
    cFLu = dot([Umean(:),-Vmean(:)]',[-fx(:),-fy(:)]',1)'.*thx(:);
    cFLv = dot([Umean(:),-Vmean(:)]',[fy(:),-fx(:)]',1)'.*thx(:);
    cFLu=reshape(cFLu,size(Umean));cFLv=reshape(cFLv,size(Umean));
    
    figure;
    imagesc(cFLu);colorbar
    
    tFL =sqrt(Umean.^2+Vmean.^2).*thx; %total flux at any point, per meter cross section!! area accoutned for with trapz below
    
%     figure;
%     imagesc(tFL);colorbar
    
%% interpolate ELA and create pixel mask
Afluxes=zeros(1,length(ELs));
fluxes=NaN;
for iEL=1:length(ELs)
M=contourc(DEM3,[ELs(iEL),ELs(iEL)]);
if numel(M)>4
[x,y]=C2xyz(M);
fluxes = zeros(1,length(x));
ELAmask= 0.*mask;

for iL=1:length(x)
    curx=x{iL}';
    cury=y{iL}';
%     curpix=[floor(cury),floor(curx);floor(cury),ceil(curx);ceil(cury),floor(curx);ceil(cury),ceil(curx)];
%     plot(curx,cury);hold on
if length(curx)>1
%     segFL=interp2(tFL,x{iL},y{iL});
    segFL=interp2(cFLu,x{iL},y{iL});
%     segTHX=interp2(thx,x{iL},y{iL});
    delx = diff(x{iL});
    dely = diff(y{iL});
    dist = sqrt(delx.^2+dely.^2);
    cdist = [0,cumsum(dist)];
%     fluxes(iL) = trapz(cdist,segFL.*segTHX);
    fluxes(iL) = trapz(cdist,segFL);
else
    fluxes(iL)=NaN;
end
%     curind=sub2ind(size(DEM),curpix(:,1),curpix(:,2));
%     ELAmask(curind)=1;
end
end
Afluxes(iEL)=nansum(fluxes);
end
% figure
% plot(ELs,Afluxes)
[FLout,iELA] = max(Afluxes);
ELA=ELs(iELA);

