function Z = imgaussfiltNaN(im,sig,varargin)
%imgaussfiltNaN - Function to implement a gaussian filter as in imgaussfilt() but ignoring NaN values. follows approach at:
%https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
%
% Syntax:  Z = imgaussfiltNaN(im,sig)
%          Z = imgaussfiltNaN(im,sig,erasenans)
%          Z = imgaussfiltNaN(______,Name,Value)
%
% Inputs:
%    im - image to be gaussian-filtered, containing NaNs
%    sig - standard deviation for gaussian kernel
%    erasenans - (optional) flag to use nan values in final output (=1, default), or to provide gap-filled filtered results (=0)
%    Name,Value - Name-Value pairs used to control aspects of the filtering as per imgaussfilt. Not tested thoroughly 
%
% Outputs:
%    Z - filtered image
%
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-Aug-2020

%% ------------- BEGIN CODE --------------
%parse inputs
nvarin=numel(varargin);
if nvarin>0
    if mod(nvarin,2)==1 %odd number of variable inputs - contains erasenans
        erasenans=varargin{1};
        if nvarin>1
            opts=[varargin{2:end}];
        end
    else % even number of inputs - erasenans not specified
        erasenans=1;
    end
else
    erasenans=1;
end



%% set up U,V inputs
nans=isnan(im); %find nans
V=im; %duplicate, replace nanas with 0
V(nans)=0;

W=1-nans; %logical array where non-nan

%% filter both
if exist('opts')
    VV=imgaussfilt(V,sig,opts);
    WW=imgaussfilt(W,sig,opts);
else
    VV=imgaussfilt(V,sig);
    WW=imgaussfilt(W,sig);
end

%ratio is filtered data ignoring NaNs
Z=VV./WW;

if erasenans==1 %option to remove NaNs again or not
    Z(nans)=NaN;
end