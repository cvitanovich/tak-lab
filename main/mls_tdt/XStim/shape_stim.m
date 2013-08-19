function [stim_out] = shape_stim(stim, ramppts, leadpts, endpts, scalefactor, ...
    itd_filt, ild_filt, hrtf_filt, ephonepts, mod_params);

%function [stim_out] = shape_stim(stim, ramppts, leadpts, endpts, scalefactor, ...
%    itd_filt, ild_filt, hrtf_filt, ephonepts, mod_params);
%
%   SHAPE_STIM: REQUIRES 5 ARGINS: stim, ramppts, leadpts, endpts, scalefactor
% and can have up to 5 additional argins: itd_filt, ild_filt, hrtf_filt, 
% ephonepts, mod_params
% mod_params should be  [mod_type mod_depth mod_freq mod_phase] or may be empty
%
% recommend scaleFactor == 6000;



if nargin < 5
    disp('SHAPE_STIM: REQUIRES 5 ARGINS: stim, ramppts, leadpts, endpts, scalefactor')
    disp('and can have up to 5 additional argins: itd_filt, ild_filt, hrtf_filt, ephonepts, mod_params');
    disp('mod_params should be  [mod_type mod_depth mod_freq mod_phase] or may be empty')
    return
end

% normalize to ACpower
stim = stim/ mom(stim,2);
% remove DC offset
stim = stim - mean(stim);
% add modulation envelope
if ~isempty('mod_params')
    DUR = length(stim)/30;
    env = make_env(DUR, mod_params(1), mod_params(2), mod_params(3), mod_params(4));
    stim = stim .* env;
end

% ramp
if ramppts
stim = [stim(1:ramppts) .* [1/ramppts:1/ramppts:1] ones(1,length(stim - 2*ramppts) ...
        stim((end-ramppts+1):end) .* [(1-1/ramppts:-1/ramppts:0]];
end

% remove DC offset
stim = stim - mean(stim);

% pad in front
if leadpts
stim = [zeros(1,leadpts) stim];
end

% pad in back
if endpts
stim = [stim zeros(1,endpts)];
end

% convolve with itd,ild,hrtf filters
if ~isempty(itd_filt);
    stim = conv(stim,itd_filt);
stim = stim - mean(stim);
end

if ~isempty(ild_filt)
    conv(stim,ild_filt);
stim = stim - mean(stim);
end

if ~isempty(hrtf_filt)
    conv(stim,ild_filt);
stim = stim - mean(stim);
end

% scale
if ~isempty('scalefactor')
stim_out = stim * scalefactor;
end

if ~isempty('ephonepts')
stim = [stim zeros(1,ephonepts)];
end    