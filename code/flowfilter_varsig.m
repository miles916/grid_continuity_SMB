function [B,misfit] = flowfilter_varsig(N,A)
%flowfilter_varsig - implements a nan-ignoring gaussian filter for variable
%distances according to 4x 'local' ice thickness
%
% Syntax:  [B,misfit] = flowfilter_varsig(N,A)
%
% Inputs:
%    N - data stack in grid_continuity_SMB structire
%    A - image to be filtered with variable distance
%
% Outputs:
%    B - filtered image
%    misfit - difference of the mean before and after filtering
%
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2020; Last revision: 20-Dec-2023

% construct special filter
fs.L = 4; %length of filter area in ice thicknesses*velocity - assume symmetrical for now
% fs.W = 1; %width of filter area in ice thicknesses
% fs.theta = atan2(N.V,N.U);


% rad1=nanmax(N.THX(N.MASK));%     find max thx, radius
% rad2=nanmean(N.Smean(N.MASK));%     find max thx, radius
% dist=3*ceil(rad1./N.DX).*ceil(rad2./N.DX);

% N.Smean2=imgaussfiltNaN(N.Smean,fs.L.*ceil(max(N.THX(:))./N.DX)); %smooth velocity signal

fs.sdist=ceil(fs.L.*N.THX./N.DX);%.*ceil(N.Smean2.*N.MASK./N.DX); %local thickness in pixel-equivalent
fs.smax=nanmax(fs.sdist(:)); %maximum distance in pixel-equivalent

% %% derive skeleton, expand, attribute nearest distance
% SK=bwskel(imbinarize(N.THX));
% %  SK = imdilate(SK,strel('diamond',10));
% SKd=fs.sdist.*SK; SKd(SKd==0)=NaN;
% [r,c]=ind2sub(size(SK),find(SK==1));
% [r2,c2]=ind2sub(size(SK),find(N.MASK==1));
% r2=1:size(SK,1);c2=1:size(SK,2);[c2,r2]=meshgrid(c2,r2);
% relD=griddata(c,r,SKd(SK),c2,r2,'natural');
% % % relD=inpaint_nans(SKd,2);
% fs.sdist=relD.*N.MASK;

% %determine maximum value in Xm radius
winLength=ceil(250./N.DX); %find within 500m window
fs.sdist = colfilt(fs.sdist, [winLength winLength], 'sliding', @max).*N.MASK;

% SE=strel('disk',8,8);%8 cell radius, for example
% nhood=SE.Neighborhood;
% f = @(x) smr(x,nhood);
% % % % test3 = colfilt(DEM,size(nhood),'sliding',f);
% testSMR=nlfilter(DEM,size(nhood),f).*dx; %SLOW

% Variable standard deviations
st = 1:fs.smax;
% Size of sliding window, Should be an odd number
m = ceil(fs.smax)+1;
kcenter = ceil(m/2);
% Image to be smoothed
% A = N.FDIV0;

% Cell with variable standard deviations of smoothing
STSIZE = int16(fs.sdist); STSIZE(STSIZE==0)=1;
st = cellfun(@(x) imgaussfiltNaN(A,x,'FilterSize',m),num2cell(st),'Unif',false);

%extract value from cells
for ix=1:numel(A)
    B(ix) = st{STSIZE(ix)}(ix);
end
B=reshape(B,size(A));
misfit= (nanmean(B(:))-nanmean(A(:)))
% B=B-(nanmean(B(:))-nanmean(A(:))); %rescale so that the integral is the same as the original

% figure;imagesc(A);caxis([-5,5])
% figure;imagesc(B);caxis([-5,5])
% N.FDIV=B;

end
