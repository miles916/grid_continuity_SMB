function N = subset_to_DH(N)

% derive new bounds
domain = abs(N.DH)>0;
domx = nansum(domain,1)>0;
domy = nansum(domain,2)>0;
ix1 = find(domx,1,'first');
ix2 = find(domx,1,'last');
iy1 = find(domy,1,'first');
iy2 = find(domy,1,'last');

%% apply bounds
N.DEM=N.DEM(iy1:iy2,ix1:ix2);
N.DH=N.DH(iy1:iy2,ix1:ix2);
N.S=N.S(iy1:iy2,ix1:ix2);
N.THX=N.THX(iy1:iy2,ix1:ix2);
N.U=N.U(iy1:iy2,ix1:ix2);
N.V=N.V(iy1:iy2,ix1:ix2);
N.x3g=N.x3g(iy1:iy2,ix1:ix2);
N.y3g=N.y3g(iy1:iy2,ix1:ix2);
N.x3=N.x3(ix1:ix2);
N.y3=N.y3(iy1:iy2);
if isfield(N,'Ue')
    N.Ue=N.Ue(iy1:iy2,ix1:ix2);
    N.Ve=N.Ve(iy1:iy2,ix1:ix2);
end
%% correct refmat
N.Rout(3,1)=N.Rout(3,1)+N.Rout(2,1).*ix1;
N.Rout(3,2)=N.Rout(3,2)+N.Rout(1,2).*iy1;