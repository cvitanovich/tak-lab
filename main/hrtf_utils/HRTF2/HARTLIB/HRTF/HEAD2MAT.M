function HEAD2MAT(file)

% HEAD2MAT merges heaphone-TF2 to matfile.
% 
%          HEAD2mat(name_core) 
%          name_core is the name of file without extension xx.sc2
%          The matlab-file will be name_core.mat

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date 26.09.1995
% Modified: 7.7.96 Klaus Hartung
%  	help improved
%       check if files with name_core exist


% Calls 	Toolbox
% fexist 	hartung/toolbox
% mtlrh		hartung/mtl
% mtlrch	hartung/mtlrch
% sizeoffile	hartung/tools


if (nargin < 1) error('Insufficient number of arguments. Name_core is missing.'); end;
if (nargin < 2) n_headphone = 1; else n_headphone=n; end;


% checking what number are available

list=[];

disp('Checking what files are available');

for n=1:99
  filename=[file  num2str(n) '.sc2'];
  if fexist(filename) list=[list n]; end;
end

if (length(list) < 1) error (['No file matches ' file ]); end;
n_headphone=length(list);
disp(['Found ' num2str(n_headphone) ' files matching ' file '.']);
disp(''); 

disp('Reading channels. This will take a while...');

filename=[file num2str(list(1)) '.sc2'] 

[comment1, filetype, info_blocks,...
           n_channels, day_of_date, month_of_date, year_of_date,...
           no_bl_per_chan, subject_no, ...
           subject_first_name, subject_last_name, ...
           n_lines, sampling_rate,...
           lo_cutoff, hi_cutoff, eff_s_to_n, no_of_aver, ... 
           no_of_qu_bits, chan_code, ampl_gain, ...
           max_ele, min_ele, status, res_type, res, ...
           channels_per_dir, map, comment2] = mtlrh(filename);

filesize=sizeoffile(filename);
n_directions = n_channels/channels_per_dir;

n_lines_v = (filesize-256*info_blocks)/(n_channels*4*2)

N_LINES=n_lines_v;

N_CHANNELS=n_headphone*2;
HEADPHONE_RAW = zeros(N_LINES,N_CHANNELS);

for n=1:n_headphone
clc
disp(['Reading channel ' num2str(n) ' of ' num2str(n_headphone) '.']);
filename=[file  int2str(list(n)) '.sc2'];


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

HPTF_RAW(:,n*2-1)=mtlrch(filename,1);
HPTF_RAW(:,n*2)  =mtlrch(filename,2);

end;

Fs = sampling_rate;
N_LINES = n_lines;
N_DIR = 0;
N_CHANNELS = n_channels; 
 
disp(['Saving data in ' file '.mat file.']);

eval(['save ' file ' HPTF_RAW N_CHANNELS N_LINES N_DIR Fs']);

% end of file head2mat.m 
