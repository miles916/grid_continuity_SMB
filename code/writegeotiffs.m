function writegeotiffs(N,THX,P)
% convert the georeferencing matrix to a map raster reference object, as required by current geotiffwrite
nrow=size(N.DEM,1);ncol=size(N.DEM,2);
Rm=N.Rout;
xlim=[Rm(3,1)+Rm(2,1)./2, Rm(3,1)+Rm(2,1)./2+Rm(2,1).*ncol]; %cell centers lie at [x y]=[row col 1]*Rm
ylim=[Rm(3,2)+Rm(1,2)./2+Rm(1,2).*nrow, Rm(3,2)+Rm(1,2)./2];
N.Rout=maprefcells(xlim,ylim,Rm(2,1),-Rm(1,2),'ColumnsStartFrom','north');
try
    Glacier=P.Glacier;
catch
    Glacier=P;
end
    % DEM
    geotiffwrite([Glacier '_DEM.tif'],uint16(N.DEM),N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % THX
    geotiffwrite([Glacier '_THX.tif'],N.THX,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % VEL
    geotiffwrite([Glacier '_Smean.tif'],N.Smean,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % dH
    geotiffwrite([Glacier '_dH.tif'],N.DH,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % FDIV
    geotiffwrite([Glacier '_FDIV.tif'],N.FDIV,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % SMB
    geotiffwrite([Glacier '_SMB.tif'],N.SMB,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % Zones
    geotiffwrite([Glacier '_zones.tif'],uint16(N.zones),N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % zFDIV
    geotiffwrite([Glacier '_zFDIV.tif'],N.z2fdiv,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    % zSMB
    geotiffwrite([Glacier '_zSMB.tif'],N.SMBz2,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))
    if P.uncertainty>0
        % FDIVu
        geotiffwrite([Glacier '_FDIVu.tif'],N.FDIVu,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

        % SMBu
        geotiffwrite([Glacier '_SMBu.tif'],N.SMBu,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

        % zFDIVu
        geotiffwrite([Glacier '_zFDIVu.tif'],N.z2fdiv_unc,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

        % zSMB
        geotiffwrite([Glacier '_zSMBu.tif'],N.SMBz2_unc,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))

    end