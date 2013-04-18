function H = earfilt(selected_frq,number_of_lines,resolution,sampling_rate)

% function H = earfilt(selected_frq,number_of_lines,resolution,sampling_rate)
%
% EARFILT evaluates a filter function of the internal ear, with a critical
%         bandwidth. The frequency response is given by the ROEX - filter from
%         Sottek R. (1993); Modelle zur Signalverarbeitung im menschlichen
%         Geh÷r; Dissertation am Lehrstuhl fuer Elektrotechnik RWTH Aachen.
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
% (p) rewritten 20.07.95 by Klaus Hartung AEA 


filter_order = 256;

% determine critical bandwidth

bandwidth = 25.0 + 75.0 * (1.0 + 1.4 * (selected_frq/1000.0)^2)^0.69;

H = [];
h = [];
b = [];

% evaluate magnitude transfer function


i=0:1:number_of_lines;
i=i/number_of_lines;
f=i*(sampling_rate/(2.0));
f_shifted=f-selected_frq;
rel_f=abs(f_shifted)/bandwidth;
prod=4.0*rel_f;

H=sqrt( (1.0 + prod) .* exp(-prod) );



b = fir2(filter_order,i,H);
sb=size(b);

h=zeros(1,1024);

h(1:sb(2))=b;

H = rfft(h);

%irplot(b',sampling_rate);
%magplot(H',sampling_rate);
%phasplot(H',sampling_rate); 
