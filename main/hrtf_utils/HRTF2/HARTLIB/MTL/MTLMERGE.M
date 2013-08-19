function mtlmerge(file_list, outputfile,uoo_flag)

% function mtlmerge(file_list, outputfile)

[n_files,str_len] = size(file_list);

filename = file_list(1,:);

[comment1, file_type_ref, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines_ref, sampling_rate_ref,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir_ref, map, comment2] = mtlrh(filename,0);



new_n_dir = 0;
total_channels = 0;

for n=1:n_files
clc
disp(['MTLMERGE: Scanning file ', num2str(n), ' of ', num2str(n_files)] );

% read header (file_list(n,:=))

filename = file_list(n,:);

[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(filename,0);

% compare header
% if different error-message and end

if (filetype ~= file_type_ref) error('different file_type'); return; end;
if (n_lines ~= n_lines_ref) error('different number_of_lines'); return; end;
if (sampling_rate ~= sampling_rate_ref) error('different sampling_rate'); return; end;
%if (first_line ~= first_line_ref) error('different first_line'); return; end;
%if (last_line ~= last_line_ref) error('different last_line'); return; end;

% compare binaural
if (channels_per_dir ~= channels_per_dir_ref) error('can t mix mono and stereo files'); return; end;

total_channels=total_channels+n_channels;

% if direction_matrix
if (exist('map')==1)
%if ((info_blocks-2) > 0)
%   map=mtlrdir(file_list(n,:));
    new_n_dir = new_n_dir +1;
end;

% set_new_direction_matrix
new_map(1:2,n)=map;

end;  % of scanning loop

% calculate header param

new_info_blocks=2+ceil((new_n_dir*2*4)/256);

channels_per_dir=total_channels/new_n_dir;


% write header(outputfile)

day_of_date= dayofdate;
mont_of_date=mofdate;
year_of_date=yofdate;


mtlwh(outputfile,...
comment1, filetype, new_info_blocks,...
total_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines_ref, sampling_rate_ref,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, new_map, comment2);



out_chan_n=1;

for n=1:n_files
filename = file_list(n,:);

[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(filename);


clc
disp(['MTLMERGE: Appending file ', num2str(n), ' of ', num2str(n_files)] );

   for c=1:n_channels
       channel=mtlrch(file_list(n,:),c);
       mtlwch(outputfile,channel,out_chan_n);
       out_chan_n = out_chan_n+1;
   end;
end;

return;

