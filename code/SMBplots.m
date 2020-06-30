function SMBplots(N,outtitle1,Glacier)
cmap = cbrewer('div','RdBu',21);
    
    cmap2 = [flipud(cbrewer('seq','Reds',11));cbrewer('seq','Blues',5)];
% cmap = [0,0,0;cbrewer('div','RdBu',21);0,0,0];

    Nout = 'grid-emergence'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(-1.*N.FDIV)
    colorbar;colormap(flipud(cmap))
    title(Ntitle)
    caxis([-10,10])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    %
    Nout = 'EL-zones'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.zones)
    colorbar;colormap([0,0,0;lines(256)])
    title(Ntitle)
%     caxis([-5,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    
    %
    Nout = 'EL-zone-avg-emergence'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(-1.*N.z2fdiv)
    colorbar;colormap(flipud(cmap))
    title(Ntitle)
    caxis([-5,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
   
    Nout = 'grid-dh'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.DH)
    colorbar;colormap(cmap2)
    title(Ntitle)
    caxis([-10,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])

    Nout = 'grid-SMB'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.SMB)
    colorbar;colormap(cmap2)
    title(Ntitle)
    caxis([-10,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    %
    Nout = 'grid-Hdensity'
    Ntitle = [Nout ' (1000 km m^{-3}), ' outtitle1];
    figure
    imagesc((N.Hdensity))
    colorbar;
    title(Ntitle)
    caxis([.5,1])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    %
    Nout = 'thx'
    Ntitle = [Nout ' (m), ' outtitle1];
    figure
    imagesc(N.THX)
    colorbar;%colormap(cmap)
    title(Ntitle)
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])

    %
    Nout = 'EL-zonal-SMB'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.SMBz2)
    colorbar;colormap(cmap2)
    title(Ntitle)
    caxis([-10,5])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])

    %
    Nout = 'surface-speed'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    figure
    imagesc(N.Smean)
    colorbar;
    title(Ntitle)
    caxis([0,100])
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
    %
    Nout = 'SMB-EL'
    Ntitle = [Nout ' (m a^{-1}), ' outtitle1];
    iN = find(N.MASK);
    figure;
    scatter(N.DEM(iN),N.SMB(iN),4,'kd');hold on
    scatter(N.DEM(iN),N.SMBz2(iN),4,'b*');
    ylabel('SMB m/a')
    ylim([-20,10])    
    title(Ntitle)
    legend('gridded SMB','zonal SMB','elevation fluxes')
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])
     %% plot SMB curves and fluxes with elevation
   
    dZ=25;
    ELS=dZ*[floor(min(N.DEM(iN)/dZ)):ceil(max(N.DEM(iN)/dZ))];
    elSMBg = zeros(1,length(ELS)-1);
    elSMBz = zeros(1,length(ELS)-1);
    elSMBe = zeros(1,length(ELS)-1);
    hyps = zeros(1,length(ELS)-1);
    for iel=1:length(ELS)-1
        cIN=N.MASK&(N.DEM>=ELS(iel))&(N.DEM<ELS(iel+1));
        hyps(iel)=nansum((cIN(:)));
        elSMBg(iel)=nanmean(N.SMB(cIN));
        elSMBz(iel)=nanmean(N.SMBz2(cIN));
        elSMBe(iel)=nanstd(N.SMB(cIN));
    end
    ELS2=unique(dZ.*floor(N.DEM(iN)./dZ));
    
    Nout = 'SMB-EL-fluxes'
    Ntitle = [Nout ', ' outtitle1];
    figure;
    yyaxis left
%     plot(ELS(1:end-1)+dZ/2,elSMBg,'b');hold on
%     plot(ELS(1:end-1)+dZ/2,elSMBz,'k');hold on
    try
        boxplot(N.SMB(iN),dZ.*floor(N.DEM(iN)./dZ),'Positions',ELS2);hold on
    catch
        boxplot(N.SMB(iN),dZ.*floor(N.DEM(iN)./dZ),'Positions',ELS2(1:end-1));hold on
    end
    ylabel('SMB m/a')
    ylim([-20,10])
    yyaxis right
    plot(N.ELs,N.ELfluxes,'linewidth',2)
    ylim([0,3*max(N.ELfluxes)])
    ylabel('Flux m^3/a')
    title(Ntitle)
%     legend('gridded','zonal','flux','location','southeast')
%     legend('flux','location','southeast')
    set(gca,'XTick',1000.*floor(min(N.DEM(iN)./1000)):500:1000.*ceil(max(N.DEM(iN)./1000)))
    set(gca,'XTickLabel',get(gca,'XTick'))
    saveas(gcf,[Glacier '_' Nout '_' outtitle1 '.png'])