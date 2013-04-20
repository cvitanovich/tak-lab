function [fg,fm]=freqgrp(fftlen,samplerate);
% FREQGRP Upper cutoff frequencies fg and center frequencies fm of
%         Zwickers frequency groups.
%         If samplerate and fftlen is given, fg and fm contain the 
%         sample number of the corresponding frequency with respect 
%         to the sample frequency and the length of the FFT (results 
%         are not exact because of limited frequency resolution!).


% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 05.07.1994 A. Raab

fg=[100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,2320,...
    2700,3150,3700,4400,5300,6400,7700,9500,12000,15500];
fm=[50,150,250,350,450,570,700,840,1000,1170,1370,1600,1850,2150,...
    2500,2900,3400,4000,4800,5800,700,8500,10500,13500];
    
if nargin>0 
  df=samplerate/fftlen; % frequency resolution
  fg=round(fg./df);     % closest sample index
  fm=round(fm./df);     
end;
 
