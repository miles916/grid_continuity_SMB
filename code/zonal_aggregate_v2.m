function [zVAR] = zonal_aggregate_v2(zones,VAR,method)
%zonal_aggregate_v2 - Function to loop through segmentation and aggregate a
%variable to a mean value ignoring NaNs
%
% Syntax:  [zVAR] = zonal_aggregate(zones,VAR)
%
% Inputs:
%    zones - segmentation (uint16 raster) of domains
%    VAR - raster of variable to be aggregated to each zone
%    method - can be 'mean' for nan-ignoring mean zonal value,'rms' for root-mean-squared (does not ignore nans) ,or 'rssn' for root-sum-squared normalized by the count (ignores nans)
%
% Outputs:
%    zVAR - mean value of VAR for each zone
%
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-Aug-2020

%% ------------- BEGIN CODE --------------
mask = isnan(zones); %find overall domain

zVAR = zeros(size(VAR));  %initialize output

%one pass over the pixels instead of one full-grid comparison per zone
inzone = ~mask;
[~,~,zidx] = unique(zones(inzone)); %zone index for every pixel in the domain
vals = double(VAR(inzone));
nz = max(zidx);
good = ~isnan(vals);

switch method
    case 'mean'
        zsum = accumarray(zidx(good),vals(good),[nz 1]);
        zcnt = accumarray(zidx(good),1,[nz 1]);
        zval = zsum./zcnt; %nan-ignoring mean
        zVAR(inzone) = zval(zidx);
    case 'rms'
        zss = accumarray(zidx,vals.^2,[nz 1]); %NaNs propagate, as rms does
        zcnt = accumarray(zidx,1,[nz 1]);
        zval = sqrt(zss./zcnt);
        zVAR(inzone) = zval(zidx);
    case 'rssn'
        zss = accumarray(zidx(good),vals(good).^2,[nz 1]);
        zcnt = accumarray(zidx(good),1,[nz 1]);
        zval = sqrt(zss)./zcnt; %root-sum-squared/N
        zVAR(inzone) = zval(zidx);
end

zVAR(mask) = NaN;