function level_difference = leveldif(l_level,r_level)

% calculates overall leveldifference  oft two spectra

% written by Klaus Hartung
% Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 
% Date 29.09.1995

l_intensity = sum( abs(l_level) )
r_intensity = sum( abs(r_level) )

level_difference =  20 * (log10(l_intensity) - log10(r_intensity));
