function filled = index_nanfill(gappy,index)
[m,n]=size(index);
gappy=double(gappy(:));
index=index(:);

%find gaps
iNA=find(isnan(gappy));

%rescale index raster to [0,100]
index2=double(index-nanmin(index));
index2=100.*index2./nanmax(index2);

%determine mean values for each index
ixv=0.5:1:99.5;
val=NaN(size(ixv));
for ix=1:length(ixv)
    val(ix)=nanmean(gappy((index2>=ixv(ix)-0.5)&(index2<=ixv(ix)+0.5)));
end

%fill in gaps with interpolation
filled=gappy;
filled(iNA)=interp1(ixv,val,index2(iNA));

filled =reshape(filled,[m,n]);
