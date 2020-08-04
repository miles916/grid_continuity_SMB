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

zonelist = unique(zones(:)); %list of zones

switch method
    case 'mean'
        for iz = 1:length(zonelist) %loop through zones and derive nanmean
            cz = zonelist(iz);
            curzone = (zones==cz);
            cval = nanmean(VAR(curzone(:)==1));
            zVAR(curzone(:)==1)=(cval);
        end
    case 'rms'
        for iz = 1:length(zonelist)
            cz = zonelist(iz);
            curzone = (zones==cz);
            cval = rms(VAR(curzone(:)==1));
            zVAR(curzone(:)==1)=(cval);
        end
    case 'rssn'
        for iz = 1:length(zonelist)
            cz = zonelist(iz);
            curzone = (zones==cz);
            select=(curzone)&(isnan(VAR)==0);
            cval = sqrt(nansum(VAR(select).^2))./sum(select(:)); %root-sum-squared/N
            zVAR(curzone(:)==1)=(cval);
        end
end

zVAR(mask) = NaN;