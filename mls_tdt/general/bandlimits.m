function [lo_limit,hi_limit] = bandlimits(center_frequency,octave_width)
%[lo_limit,hi_limit] = bandlimits(center_frequency,octave_width)
%Returns the bandpass limits for a given center frequency and a given width in terms of octave
%Center frequency is in Hz

lo_limit = 2 * center_frequency / (1 + 2^octave_width);
hi_limit = lo_limit * 2^octave_width;

return