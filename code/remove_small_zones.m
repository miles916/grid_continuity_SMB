function zonesout=remove_small_zones(zonesin,mask,N)
%remove_small_zones - Function to remove small, isolated zones from a
%segmented glacier mask
%
% Syntax:  zonesout=remove_small_zones(zonesin,mask,N)
%
% Inputs:
%    zonesin - preceding segmentation (uint16 raster) including small zones
%    mask - binary raster of glacier area
%    N - minimum number of pixels for a given zone; gets overwritten for
%    small glaciers
%
% Outputs:
%    zonesout - uint16 raster of glacier segmentation with small zones eliminated
%
% Other m-files required: inpaint_nans.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-June-2020

%% ------------- BEGIN CODE --------------
% Ls=unique(zonesin(:));

if sum(mask(:))<500 %if the glacier is small, overwrite N
    N=10;
end
if sum(mask(:))>50 %if very small, skip altogether
    z1=zonesin.*uint16(mask); %ensure cropped to mask

    %find labels of current zones
    ulabs = unique(zonesin(:)); 
    ulabs=ulabs(ulabs>0);

    %% find size of each zone 
    sizes= double(z1.*0); %initialize
    for iz = 1:length(ulabs)
        curdata = z1==ulabs(iz); %find pixels in current zone
        curlabs = bwlabel(curdata); 
        curstats = regionprops(curdata,'Area');
        cursizes=0*curdata;
        curAs = [curstats.Area];
        if max(curlabs(:)>0)
            cursizes(curlabs>0) = curAs(curlabs(curlabs>0)); %relabel current zone based on the number of pixels
        end
        sizes = sizes+cursizes;
    end


    %remove small zones
    z2=z1; %initialize modified zones
    z2(sizes<N)=0; %remove those with small areas

    [~,idx] = bwdist(z2>0); %find nearest zone (euclidean distance) for all unlabeled zones

    z3=z2;
    z3(z3==0)=z3(idx(z3==0)); %relabel erased zones with nearest zone

    %crop to mask
    zonesout = uint16(z3).*uint16(mask); 
else
    zonesout=zonesin;
end
