function filled = index_nanfill(gappy,index)
[m,n]=size(index);
gappy=double(gappy(:));
index=index(:);

%find gaps
iNA=find(isnan(gappy)|gappy==0);

%rescale index raster to [0,100]
index2=double(index-nanmin(index));
index2=100.*index2./nanmax(index2);

%determine mean values for each index
ixv=0.5:1:99.5;
val=NaN(size(ixv));
for ix=1:length(ixv)
    cur=gappy((index2>=ixv(ix)-0.5)&(index2<=ixv(ix)+0.5));
    val(ix)=nanmedian(cur((cur==0)==0));
end

irem=isnan(val);
ixv(irem)=[];
val(irem)=[];


%fill in gaps with interpolation
filled=gappy;
filled(iNA)=interp1(ixv,val,index2(iNA),'linear','extrap');

filled =reshape(filled,[m,n]);
