function Dvec = movingslope(vec,supportlength,varargin)
% adapted from see https://ch.mathworks.com/matlabcentral/fileexchange/16997-movingslope
% movingslope: estimate local slope for a sequence of points, using a
% sliding window,  but ignoring NaN values. follows approach at:
%https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
% usage: Dvec = movingslope(vec)
% usage: Dvec = movingslope(vec,supportlength)
% usage: Dvec = movingslope(vec,supportlength,modelorder)
% usage: Dvec = movingslope(vec,supportlength,modelorder,dt)
%
%
% movingslope uses filter to determine the slope of a curve stored
% as an equally (unit) spaced sequence of points. A patch is applied
% at each end where filter will have problems. A non-unit spacing
% can be supplied.
%
% Note that with a 3 point window and equally spaced data sequence,
% this code should be similar to gradient. However, with wider
% windows this tool will be more robust to noisy data sequences.
%
%
% arguments: (input)
%  vec - row of column vector, to be differentiated. vec must be of
%        length at least 2.
%
%  supportlength - (OPTIONAL) scalar integer - defines the number of
%        points used for the moving window. supportlength may be no
%        more than the length of vec.
%
%        supportlength must be at least 2, but no more than length(vec)
%
%        If supportlength is an odd number, then the sliding window
%        will be central. If it is an even number, then the window
%        will be slid backwards by one element. Thus a 2 point window
%        will result in a backwards differences used, except at the
%        very first point, where a forward difference will be used.
%
%        DEFAULT: supportlength = 3
%
%  modelorder - (OPTIONAL) - scalar - Defines the order of the windowed
%        model used to estimate the slope. When model order is 1, the
%        model is a linear one. If modelorder is less than supportlength-1.
%        then the sliding window will be a regression one. If modelorder
%        is equal to supportlength-1, then the window will result in a
%        sliding Lagrange interpolant.
%
%        modelorder must be at least 1, but not exceeding
%        min(10,supportlength-1)
%
%        DEFAULT: modelorder = 1
%
%  dt - (OPTIONAL) - scalar - spacing for sequences which do not have
%        a unit spacing.
%
%        DEFAULT: dt = 1
%
% arguments: (output)
%  Dvec = vector of derivative estimates, Dvec will be of the same size
%        and shape as is vec.
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
    VV=movingslope(V,opts);
    WW=movingslope(W,opts);
else
    VV=movingslope(V);
    WW=movingslope(W);
end

%ratio is filtered data ignoring NaNs
Z=VV./WW;

if erasenans==1 %option to remove NaNs again or not
    Z(nans)=NaN;
end