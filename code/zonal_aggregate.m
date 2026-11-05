function [zVAR] = zonal_aggregate(zones,VAR)
%zonal_aggregate - Function to loop through segmentation and aggregate a
%variable to a mean value ignoring NaNs
%
% Syntax:  [zVAR] = zonal_aggregate(zones,VAR)
%
% Inputs:
%    zones - segmentation (uint16 raster) of domains
%    VAR - raster of variable to be aggregated to each zone
%
% Outputs:
%    zVAR - mean value of VAR for each zone
%
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-June-2020

%% ------------- BEGIN CODE --------------
mask = isnan(zones); %find overall domain

zVAR = zeros(size(VAR)); %initialize output

%one pass over the pixels instead of one full-grid comparison per zone
inzone = ~mask;
[~,~,zidx] = unique(zones(inzone)); %zone index for every pixel in the domain
vals = double(VAR(inzone));
good = ~isnan(vals);
nz = max(zidx);
zsum = accumarray(zidx(good),vals(good),[nz 1]); %sum omitting NaNs
zcnt = accumarray(zidx(good),1,[nz 1]); %count omitting NaNs
zmean = zsum./zcnt; %NaN for zones with no valid values, as nanmean

zVAR(inzone) = zmean(zidx); %write into new file

zVAR(mask) = NaN; %erase 0s outside mask