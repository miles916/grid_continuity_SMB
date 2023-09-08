function writegeotiffs(N,THX,P)
if THX.info.GeoTIFFTags.GeoKeyDirectoryTag.GTRasterTypeGeoKey==1 %Matlab geotiffexport off by 1 pixel
    N.Rout(3,1)=N.Rout(3,1)-N.DX;
    N.Rout(3,2)=N.Rout(3,2)+N.DX;
end
try
    Glacier=P;
catch
    Glacier=P.Glacier;
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