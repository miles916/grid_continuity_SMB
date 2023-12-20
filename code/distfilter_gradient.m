function [Bx,By] = distfilter_gradient(N,A)
%distfilter_gradient - implements a variable-range determination of x and y
%gradients overdistances according to 4x 'local' ice thickness
%
% Syntax:  [Bx,By] = distfilter_gradient(N,A)
%
% Inputs:
%    N - data stack in grid_continuity_SMB structire
%    A - image for gradients with variable distance
%
% Outputs:
%    Bx - gradient of A in x direction
%    By - gradient of A in y direction
%
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2023; Last revision: 20-Dec-2023

% construct special filter
fs.L = 4; %length of filter area in ice thicknesses*velocity - assume symmetrical for now
% fs.W = 1; %width of filter area in ice thicknesses
% fs.theta = atan2(N.V,N.U);

fs.sdist=ceil(fs.L.*N.THX./N.DX);%.*ceil(N.Smean2.*N.MASK./N.DX); %local thickness in pixel-equivalent
fs.smax=nanmax(fs.sdist(:)); %maximum distance in pixel-equivalent

% %determine maximum value in Xm radius
winLength=ceil(250./N.DX); %find within 500m window
fs.sdist = colfilt(fs.sdist, [winLength winLength], 'sliding', @max).*N.MASK;

% Variable ranges
dAdy = 2:fs.smax;dAdx=dAdy;
% Size of sliding window, Should be an odd number
m = ceil(fs.smax)+1;
kcenter = ceil(m/2);
% Image to be smoothed
% A = N.FDIV0;

% Cell with variable standard deviations of smoothing
STSIZE = int16(fs.sdist); STSIZE(STSIZE<=1)=2;%minimum distance for sliding filter

A(isnan(A))=0; % otherwise invalidates line processing

Ay=reshape(A,[1,numel(A)]);
dAdy = cellfun(@(x) movingslopeNaN(Ay,x,1,-N.DX),num2cell(dAdy),'Unif',false); %negative because y axis is flipped

Ax=reshape(A',[1,numel(A)]);
dAdx = cellfun(@(x) movingslopeNaN(Ax,x,1,N.DX),num2cell(dAdx),'Unif',false); 

%extract value from cells
for ix=1:numel(Ay)
    Bx(ix) = dAdx{STSIZE(ix)-1}(ix);%-1 is because we start at 2
    By(ix) = dAdy{STSIZE(ix)-1}(ix);
end
Bx=reshape(Bx,size(A'))';
By=reshape(By,size(A));


figure;imagesc(Bx.*N.MASK);colorbar
figure;imagesc(By.*N.MASK);colorbar

end
