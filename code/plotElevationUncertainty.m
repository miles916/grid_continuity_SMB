function plotElevationUncertainty(N,Glacier,outtitle1)
if isfield(N,'ELEVBINS')
    binGrid = N.ELEVBINS;
else
    binGrid = 50 .* floor(N.DEM ./ 50);
end
binValues = unique(binGrid(N.MASK & isfinite(binGrid)));
if isempty(binValues)
    return
end

smbMean = NaN(size(binValues));
smbUncertainty = NaN(size(binValues));
fdivMean = NaN(size(binValues));
fdivUncertainty = NaN(size(binValues));
for iBin = 1:numel(binValues)
    inBin = N.MASK & binGrid == binValues(iBin);
    smbMean(iBin) = mean(N.SMB(inBin),'omitnan');
    fdivMean(iBin) = mean(N.FDIV(inBin),'omitnan');
    nSmb = nnz(inBin & isfinite(N.SMBu));
    nFdiv = nnz(inBin & isfinite(N.FDIVu));
    if nSmb > 0
        if isfield(N,'resultMode') && strcmp(N.resultMode,'elevation_bins')
            smbUncertainty(iBin) = mean(N.SMBu(inBin),'omitnan');
        else
            smbUncertainty(iBin) = sqrt(sum(N.SMBu(inBin).^2,'omitnan')) ./ nSmb;
        end
    end
    if nFdiv > 0
        if isfield(N,'resultMode') && strcmp(N.resultMode,'elevation_bins')
            fdivUncertainty(iBin) = mean(N.FDIVu(inBin),'omitnan');
        else
            fdivUncertainty(iBin) = sqrt(sum(N.FDIVu(inBin).^2,'omitnan')) ./ nFdiv;
        end
    end
end

figure
subplot(2,1,1)
errorbar(binValues,smbMean,smbUncertainty,'o-','LineWidth',1)
ylabel('SMB (m i.e. yr^{-1})')
grid on
title(['Elevation-bin SMB with Monte Carlo uncertainty, ' outtitle1])
subplot(2,1,2)
errorbar(binValues,fdivMean,fdivUncertainty,'o-','LineWidth',1)
xlabel('Elevation-bin left edge (m a.s.l.)')
ylabel('FDIV (m i.e. yr^{-1})')
grid on
title('Elevation-bin FDIV with Monte Carlo uncertainty')
saveas(gcf,[Glacier '_SMB-FDIV-uncertainty-EL_' outtitle1 '.png'])
end
