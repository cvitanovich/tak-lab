function sc12mat(file)

% sc12mat converts sc1 file to matfile
% mtl2mat(file) file is name of mtlfile without extension .mtl
% The matlab-file will be file.mat

% written by Klaus Hartung
% Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 
% Date 26.09.1995


disp('Reading channels. This will take a while...');

filename=[file '1.sc1'];


%[filetype, info_blocks,...
%n_channels, n_lines, sampling_rate,...
%first_line, last_line, n_directions,...
%comment1, comment2] = mtlrh(filename);

           [comment1, filetype, info_blocks,...
           n_channels, day_of_date, month_of_date, year_of_date,...
           no_bl_per_chan, subject_no, ...
           subject_first_name, subject_last_name, ...
           n_lines, sampling_rate,...
           lo_cutoff, hi_cutoff, eff_s_to_n, no_of_aver, ... 
           no_of_qu_bits, chan_code, ampl_gain, ...
           max_ele, min_ele, status, res_type, res, ...
           channels_per_dir, map, comment2] = mtlrh(filename)

filesize=sizeoffile(filename);

n_lines_v = (filesize-256*info_blocks)/(n_channels*4*2);
N_LINES=n_lines;

N_CHANNELS=11;
HRTF_RAW = zeros(N_LINES/2+1,N_CHANNELS);
Fs = sampling_rate;
N_LINES = n_lines;
N_DIR = 11;
N_CHANNELS = n_channels; 

for n=1:N_DIR
 % filename=[file '_' int2str(n) '.sc2'];
  filename=[file int2str(n) '.sc1'];


  % load header information


           [comment1, filetype, info_blocks,...
           n_channels, day_of_date, month_of_date, year_of_date,...
           no_bl_per_chan, subject_no, ...
           subject_first_name, subject_last_name, ...
           n_lines, sampling_rate,...
           lo_cutoff, hi_cutoff, eff_s_to_n, no_of_aver, ... 
           no_of_qu_bits, chan_code, ampl_gain, ...
           max_ele, min_ele, status, res_type, res, ...
           channels_per_dir, map, comment2] = mtlrh(filename);
%future extension
%MAP(:,n)=map

LSTF_RAW(:,n)=mtlrch(filename,1);

end;


clear  filetype info_blocks n_channels lo_cutoff hi_cutoff  n_directions sampling_rate;


clear filename n 
filename='LS';
disp('Saving data in *.mat file.');

eval(['save ' filename ' Fs N_LINES N_DIR N_CHANNELS LSTF_RAW']);
 
