function [ELA,FLout,cFLu,tFL,ELs,Afluxes]=through_fluxes(mask,DEM,Umean,Vmean,thx,dx)%,dy,sig_H,UE,VE,ERRORs)

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

%% calculate flux corrected for slope direction (at ELA
%     FLx=Umean.*thx;
%     FLy=Vmean.*thx;
    
    cFLu = dot([Umean(:),-Vmean(:)]',[-fx(:),-fy(:)]',1);
    cFLv = dot([Umean(:),-Vmean(:)]',[fy(:),-fx(:)]',1);
    cFLu=reshape(cFLu,size(Umean));cFLv=reshape(cFLv,size(Umean));
    
%     figure;
%     imagesc(cFLu.*thx);colorbar
    
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
    segFL=interp2(tFL,x{iL},y{iL});
%     segFL=interp2(cFLu,x{iL},y{iL});
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

