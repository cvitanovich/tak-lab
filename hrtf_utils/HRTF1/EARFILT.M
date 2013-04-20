function H = earfilt(selected_frq,number_of_lines,resolution,sampling_rate)

% function H = earfilt(selected_frq,number_of_lines,resolution,sampling_rate)
%
% EARFILT evaluates a filter function of the internal ear, with a critical
%         bandwidth. The frequency response is given by the ROEX - filter from
%         Sottek R. (1993); Modelle zur Signalverarbeitung im menschlichen
%         Gehör; Dissertation am Lehrstuhl fuer Elektrotechnik RWTH Aachen.
%         With the frequency response a 100 th order FIR digital filter is
%         designed.
%
%         meaning of the input arguments:
%         - selected_frq: middle frequency of the bandpass
%         - number_of_lines: number of lines for which the filter must be
%                            designed
%         - resolution: resolution of the frequency lines
%         - sampling_rate
%         meaning of the output argument:
%         - H: complex transfer function of the internal ear with the
%              bandpass middle frequency 'selected_frq'
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter


% determine critical bandwidth

bandwidth = 25 + 75 * (1 + 1.4 * (selected_frq/1000)^2)^0.69;

% evaluate magnitude range

H = [];
for i=1:number_of_lines
   f = (i-1) * resolution;
   f = f - selected_frq;              % Verschiebung der Mittenfrequenz
   H(i) = sqrt( (1 + 4*abs(f)/bandwidth)*exp(-4*abs(f)/bandwidth) );
end;

% evaluate phase range

filter_order = 100;
frq_range = linspace(0,sampling_rate/2,number_of_lines);
frq_range = frq_range / (sampling_rate/2);

b = [];
b = fir2(filter_order,frq_range,H);
b = [b,zeros(1,2*(number_of_lines-1) - filter_order - 1)];

H = fft(b);
H = H(1:number_of_lines);

% plotten der Übertragungsfunktion des Filters
%figure;
%subplot(2,1,1);
%h = real(H);
%plot(h);
%subplot(2,1,2);
%ph = imag(H);
%plot(ph);
%pause;ÿ
