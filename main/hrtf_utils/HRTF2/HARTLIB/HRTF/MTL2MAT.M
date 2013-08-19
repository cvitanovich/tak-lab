function mtl2mat(file)

% MTL2MAT converts mtl2 file to matfile
% mtl2mat(file) file is name of mtlfile without extension .mtl
% The matlab-file will be file.mat

% written by Klaus Hartung
% Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 
% Date 26.09.1995


filename=[file '.mtl'];
%filename=[file ]
disp('Open .mtl file.');

% load header information

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
           channels_per_dir, MAP, comment2] = mtlrh(filename);

filesize=sizeoffile(filename);

%n_lines_v = (filesize-256*info_blocks)/(n_channels*4*2)
%n_lines=n_lines_v;

Fs = sampling_rate;
N_LINES = n_lines;
N_DIR = n_channels/channels_per_dir;
N_CHANNELS = n_channels; 

clear   info_blocks n_channels first_line last_line n_directions sampling_rate;

% load directions 

%MAP = mtlrdir(filename)

% load channels
if (filetype==3)
  HRTF_RAW = zeros(N_LINES/2+1,N_CHANNELS);
end
if (filetype==2)
  HRIR_RAW= zeros(N_LINES,N_CHANNELS);
end

disp('Reading channels. This will take a while...');

for n=1:N_CHANNELS
  HRTF_RAW(:,n)=mtlrch(filename,n);
end

% save as matlab file
disp('Saving File');

clear filename n 

disp('Saving data in *.mat file.');

eval(['save ' file ' Fs N_LINES N_DIR N_CHANNELS MAP HRTF_RAW ']);
 
