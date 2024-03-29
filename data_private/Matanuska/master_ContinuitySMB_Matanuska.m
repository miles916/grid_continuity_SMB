%master_ContinuitySMB_ex.m - Master script to estimate surface mass balance
%distribution for a glacier from inputs of ice thickness, thinning, and
%velocity, based on the continuity equation (see, e.g. Miles et al, 2021)
%
% New glaciers require this template to be coded with paths for the full set of
% inputs, and any additional needed preprocessing. 
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-June-2020

%% ------------- BEGIN CODE --------------
clear 
close all

%set home directory
homedir = 'C:\Users\miles\Documents\GitHub\grid_continuity_SMB';

%glacier ID
P.Glacier = 'Matanuska'

%title of output directory
datatitle = ['test_' P.Glacier];

%inputs        
V.pathx = fullfile(homedir,'data_private',P.Glacier,'Vx_Matanuska_ITSLIVE.tif');
V.pathy = fullfile(homedir,'data_private',P.Glacier,'Vy_Matanuska_ITSLIVE.tif');
V.mult=1; %scale to convert input units to m/a
P.Vreproj=1; %0 if velocities are provided oriented in the correct coordinate system, 1 if they are in the source data projection, [2 if they are true north], [3 to determine from slope]
DH.path = fullfile(homedir,'data_private',P.Glacier,'dhdt_Matanuska_2007_2016.tif');
DH.mult=1; %scale to convert input units to m/a
THX.path = fullfile(homedir,'data_private',P.Glacier,'RGI60-01.10557_thickness_composite.tif');
DEM.path = fullfile(homedir,'data_private',P.Glacier,'GDEM_Matanuska.tif');

%settings
P.plotouts=1; %output plots or not
P.exports=1; %save geotiffs or not
P.DX = 250; %resolution to run calculations at 
P.segdist=500; %effective linear distance between flowbands
P.Vfilter=0; %swtich to smooth velocity data or not
P.dhfilter=0; %switch to peform 2x 3-sigma outlier removal (from overall dataset - only if erroneous pixels are common)
P.THXfilter=0; %switch to apply simply Gaussian filter to thickness data (useful for field thickness measurements)
P.umult=0; %switch for column-average velocity [0.8-1] is physical range. '0' estimates based on THX dist. '2' estimates for each pixel. 
% P.Vreproj=0; %0 if velocities are provided oriented in the correct coordinate system, 1 if they are in the source data projection, [2 if they are true north], [3 to determine from slope]
P.fdivfilt=2; %use VanTricht gradient filters (2), just flux filter (1) or not at al (0)
P.uncertainty=0; %0 for 'simple run without uncertainty, otherwise N runs; note that you will then need to setup perturbations within the N structure

% initialization
addpath(genpath([homedir '\code'])) %add path to related scripts

outdir = [homedir '\results\' datatitle '_' date]
mkdir(outdir)

outtitle1 = [num2str(P.DX) 'mgrid'];
cd(outdir)
    
    %% Load all input data, resample to common grid, etc
    [THX,DEM,V,DH]=load_subset_all(P,THX,DEM,V,DH);
    
    %% Cleaning raw data if needed
    
    %identify likely errors
    V.iERR=(abs(V.Uraw)>400)|(abs(V.Vraw)>400);
    if P.Vfilter==1 %gaussian low-pass filter removing extreme variations
        V.Uraw=imgaussfilt(V.Uraw,5);
        V.Vraw=imgaussfilt(V.Vraw,5);
    end

    if P.THXfilter==1% filter thickness data
        THX.dataR=THX.data;
        THX.data=imgaussfilt(THX.dataR,4); %gaussian low-pass filter. Important for thickness maps derived from field data
        THX.data(MASK0==0)=0;
    end
    
    if P.dhfilter==1 %if noisy or gappy, apply some preprocessing to the dH
        DH.data2=DH.data;
        DH.errthresh1=3.*nanstd(DH.data(:));
        DH.data2(abs(DH.data2)>DH.errthresh1)=NaN;
        DH.errthresh2=3.*nanstd(DH.data2(:));
        DH.data2(abs(DH.data2)>DH.errthresh2)=NaN;
        DH.dH3=imgaussfilt(DH.data2); %smooths dH slightly, expands NaN around bad DH data
    else
        DH.dH3=DH.data;
    end
    
    %% REPROJECT AND RESAMPLE ALL INPUTS (to THX coordinate system), derive MASK, etc
    N = resample_inputs(P.DX,THX,V,DEM,DH);
    N.fdivfilt=P.fdivfilt;
    N.uncertainty=P.uncertainty;
    
    %% fill gaps in reprojected inputs if needed
    
    %gap-fill dH based on elevation
    N.DH0=N.DH; %unfilled values
    N.DH(N.MASK)= index_nanfill(N.DH(N.MASK),N.DEM(N.MASK));
    N.DH((N.MASK==0))=0;

    %% zone segmentation
%     N.zones = uint16(segment_Gmask_slope2(N.DEM,N.MASK,P.DX,P.segdist)); %segments glacier mask and DEM into approximately uniformly-spaced elevation bands
    N.zones = uint16(segment_Gmask_EL(N.DEM,N.MASK,P.DX,P.segdist./5)); %segments glacier mask and DEM into approximately uniformly-spaced elevation bands

    %% determine column-average velocity correction
    N = derive_umult(N,P);

    %% SMB calculations
    cd(outdir)

    if N.uncertainty>1
        N=FluxCalcsUncertainty(N); %calculates fluxes 
    else
        N=FluxCalcsSimple(N); %calculates fluxes 
    end    
%% DETERMINE FLUXES THROUGH EACH ELEVATION BAND
%     [N.ELA,N.FLout,N.cFLu,N.tFL,N.ELs,N.ELfluxes]=through_fluxes(N.MASK,N.DEM,N.Umean,N.Vmean,N.THX,N.DX);%,dy,sig_H,UE,VE,ERRORs)

if P.plotouts==1
    SMBplots(N,outtitle1,P.Glacier);
end

    
    %% export geotifs
if P.exports==1
    writegeotiffs(N,THX,P)
end

%%    

disp('finished')
   
