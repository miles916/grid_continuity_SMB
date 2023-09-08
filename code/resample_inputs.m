function N = resample_inputs(DX,THX,V,DEM,DH)
    buffdist = 1000; %expands the domain for subsetting around the glacier

    %new coordinates based on the THX data, but at DX interval
    N.x3 = DX.*[(floor((THX.xm(1)-buffdist)/DX)):((ceil(THX.xm(end)+buffdist)/DX))];
    N.y3 = DX.*[(ceil((THX.ym(1)+buffdist)/DX)):-1:(floor((THX.ym(end)-buffdist)/DX))];
    [N.x3g,N.y3g] = meshgrid(N.x3,N.y3);

    N.Rout = [0,-DX;DX,0;N.x3(1)-DX,N.y3(1)+DX]; %updated georeferencing matrix based on new domain
    N.DX=DX;

    %RESAMPLE THICKNESS
    THX.data(isnan(THX.data))=0;%Nans result in contraction of the glacier domain after interpolation
%     [THX.xN,THX.yN] = projfwd(THX.info,THX.LatG,THX.LonG); %compute projected coordinates
%     N.THX = griddata(THX.xN(:),THX.yN(:),double(THX.data(:)),N.x3g(:),N.y3g(:),'cubic');
    N.THX = griddata(THX.xmG(:),THX.ymG(:),double(THX.data(:)),N.x3g(:),N.y3g(:),'linear');
    N.THX = reshape(N.THX,size(N.x3g));
    N.MASK=N.THX>0;
    
    %resample velocity data
    [V.xN,V.yN] = projfwd(THX.info,V.LatG(:),V.LonG(:)); %compute projected coordinates
    V.xN(V.iERR)=[];V.yN(V.iERR)=[];
    if V.Vreproj==0 %if no reprojection, assume north-oriented
        V.Uraw2=V.Uraw(:);V.Uraw2(V.iERR)=[]; %OLD
        V.Vraw2=V.Vraw(:);V.Vraw2(V.iERR)=[]; %OLD
    else %use projected end-points
        [V.x2N,V.y2N] = projfwd(THX.info,V.Lat2(:),V.Lon2(:)); %compute projected coordinates for vector end-points
        V.x2N(V.iERR)=[];V.y2N(V.iERR)=[];
        V.Uraw2=V.x2N-V.xN; %velocity vector in new coord system
        V.Vraw2=V.y2N-V.yN; %velocity vector in new coord system
    end
    N.U = griddata(V.xN(:),V.yN(:),double(V.Uraw2(:)),N.x3g(:),N.y3g(:),'linear');
    N.U = reshape(N.U,size(N.x3g));
    N.V = griddata(V.xN(:),V.yN(:),double(V.Vraw2(:)),N.x3g(:),N.y3g(:),'linear');
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
