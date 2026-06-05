% TSinterp = interp_over_nan(T,TS,method)
% 
% Simple linear interpolation over a time series (or similar data set) TS
% with a corresponding grid T. The purpose of this is to eliminate NaN's in
% favor of a linear interpolation, something linterp doesn't do by default.
% If there's a NaN at the beginning and/or end, this script simply makes it
% the same as the nearest non-Nan number.
% 
% method is a string, one of:
% ['nearest', 'next', 'previous', 'linear', 'spline', 'pchip', 'cubic']
% and is the method used by MATLAB's built-in interp1 function, called in
% this function. Deafault is 'linear', just like interp1.


function TSinterp = interp_over_nan(T,TS,varargin)

if nargin == 2
    method = 'linear';
elseif nargin == 3
    method = varargin{1};
else
    error('Only 2 or 3 arguments are expected.')
end

TSnonan = TS(isfinite(TS));
% get rid of leading and trailing NaN's:
if isempty(TSnonan)
    TSinterp = TS;
else
    TS(1) = TSnonan(1);
    TS(end) = TSnonan(end);

    TSnonan = TS(isfinite(TS));
    Tnonan = T(isfinite(TS));

    TSinterp = interp1(Tnonan,TSnonan,T,method);
end

end
