function Nout = upsample_results(N, info)
%UPSAMPLE_RESULTS Expand coarsened result rasters back to the native grid.
% Nearest-neighbour (block replication), masked to the native glacier mask.

Nout = N;
Nout.DX = info.nativeDX;
if info.factor <= 1
    return
end
f = info.factor;
sz = info.nativeSize;
fields = {'SMB', 'FDIV', 'Hdensity', 'SMBu', 'FDIVu'};
for i = 1:numel(fields)
    if isfield(N, fields{i})
        A = repelem(N.(fields{i}), f, f);
        A = A(1:sz(1), 1:sz(2));
        A(~info.nativeMask) = NaN;
        Nout.(fields{i}) = A;
    end
end
Nout.MASK = info.nativeMask;
end
