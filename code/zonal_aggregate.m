function [zVAR] = zonal_aggregate(zones,VAR)

mask = isnan(zones);

zVAR = zeros(size(VAR));

zonelist = unique(zones(:));

for iz = 1:length(zonelist)
	cz = zonelist(iz);
	curzone = (zones==cz);
	cfdiv = nanmean(VAR(curzone(:)==1));
	zVAR(curzone(:)==1)=(cfdiv);
end

zVAR(mask) = NaN;