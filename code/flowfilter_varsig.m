function [B,misfit] = flowfilter_varsig(N,A)

% construct special filter
fs.L = 4; %length of filter area in ice thicknesses*velocity - assume symmetrical for now
% fs.W = 1; %width of filter area in ice thicknesses
% fs.theta = atan2(N.V,N.U);


% rad1=nanmax(N.THX(N.MASK));%     find max thx, radius
% rad2=nanmean(N.Smean(N.MASK));%     find max thx, radius
% dist=3*ceil(rad1./N.DX).*ceil(rad2./N.DX);

% N.Smean2=imgaussfiltNaN(N.Smean,fs.L.*ceil(max(N.THX(:))./N.DX)); %smooth velocity signal

fs.sdist=fs.L.*ceil(N.THX./N.DX);%.*ceil(N.Smean2.*N.MASK./N.DX); %local thickness in pixel-equivalent
fs.smax=nanmax(fs.sdist(:)); %maximum distance in pixel-equivalent

% SE=strel('disk',8,8);%8 cell radius, for example
% nhood=SE.Neighborhood;
% f = @(x) smr(x,nhood);
% % % % test3 = colfilt(DEM,size(nhood),'sliding',f);
% testSMR=nlfilter(DEM,size(nhood),f).*dx; %SLOW

% Variable standard deviations
st = 1:fs.smax;
% Size of sliding window, Should be an odd number
m = fs.L*ceil(fs.smax)+1;
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
