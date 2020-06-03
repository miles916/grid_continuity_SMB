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

zonelist = unique(zones(:)); %list of zones

%loop through zones and derive nanmean
for iz = 1:length(zonelist)
	cz = zonelist(iz); %identify current zone
	curzone = (zones==cz); %current domain
	cVAR = nanmean(VAR(curzone(:)==1)); %mean value omitting NaNs
	zVAR(curzone(:)==1)=(cVAR); %write into new file
end

zVAR(mask) = NaN; %erase 0s outside mask