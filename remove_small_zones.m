function zonesout=remove_small_zones(zonesin,mask,N)
%determine size of zones (anything less than 50 pix is small)
% Ls=unique(zonesin(:));

if sum(mask(:))<500
    N=10;
end
if sum(mask(:))>50
    z1=zonesin.*uint16(mask);

    ulabs = unique(zonesin(:));
    ulabs=ulabs(ulabs>0);

    sizes= double(z1.*0);
    for iz = 1:length(ulabs)
        curdata = z1==ulabs(iz); %find current elements
        curlabs = bwlabel(curdata);
        curstats = regionprops(curdata,'Area');
        cursizes=0*curdata;
        curAs = [curstats.Area];
        if max(curlabs(:)>0)
            cursizes(curlabs>0) = curAs(curlabs(curlabs>0));
        end
        sizes = sizes+cursizes;
    end


    %remove small zones
    % z2=double(z1);
    z2=z1;
    z2(sizes<N)=0;

    % removeIDX = find([stats.Area]<N);
    % z2(ismember(z2(:)',removeIDX))=0;

    %interpolate with nearest neighbors
    % z3 =fillmissing(z2,'nearest'); %works only in columns or rows

    % z2(mask==0)=NaN;
    % z3= inpaint_nans(z2); %gives decimals...
    [~,idx] = bwdist(z2>0);

    z3=z2;
    z3(z3==0)=z3(idx(z3==0));

    %crop to mask
    zonesout = uint16(z3).*uint16(mask);
else
    zonesout=zonesin;
end
