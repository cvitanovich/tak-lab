function level_difference = leveldif(left_X,right_X)

% function level_difference = leveldif(left_X,right_X)
%
% LEVELDIF evaluates the interaural level difference by the following method.
%          The transfer functions of the external ears are filtered by a 
%          ROEX-filter (see EARFILT) with the critical bandwidth depending on
%          the choosen frequency to filter. The intensities of the the left and
%          the right channel are evaluated and the level difference related
%          to the left ear is evaluated and returned in [dB].
%
%          Description of the input arguments:
%          - left_X: the filtered transfer function of the left ear
%          - right_X: the filtered transfer function of the right ear
%
%          For further information see: ia_eval1, timediff, earfilt
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


l_intensity = sum( (abs(left_X)).^2 );
r_intensity = sum( (abs(right_X)).^2 );

level_difference =  10 * log10(l_intensity) - 10 * log10(r_intensity);