function mtlplot(catalogue, elevation, azimuth, lp)

% mtlplot(catalogue, elevation, azimuth,lp)
%
% mtlplot(catalogue, channel,lp)
%
% mtlplot(catalogue, lp)
%
% mtlplot(catalogue)

% MTLPLOT plots the impulse response or transfer function of the selected
% direction or channel

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995

min_frq=200;
max_frq=15000;

[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, direction_matrix, comment2] = mtlrh(catalogue);

%direction_matrix=mtlrdir(catalogue);





if (nargin == 1) 
	for index=1:n_channels
		first_line=min_frq*(n_lines/sampling_rate);
		last_line=max_frq*(n_lines/sampling_rate);
		channel=mtlrch(catalogue,index);
                if (filetype==3)
	                l=abs2db(channel(first_line:last_line));
                 	max_level(index)=max(l);
			min_level(index)=min(l);
                end;
	end;

	max_abs=max(max_level);
	min_abs=min(min_level);


	for index=1:n_directions

           elevation=direction_matrix(1,index);
           azimuth=direction_matrix(2,index);
           
           title=sprintf('Catalogue : %s, Elevation: %i , Azimuth: %i', catalogue, elevation,azimuth);

	   channel=mtlrch(catalogue,index);
 
   	   if (channels_per_dir==2)
                g=10^(max_abs/20);
        	channel1=mtlrch(catalogue,index*2-1)/g;
	        channel2=mtlrch(catalogue,index*2)/g;
	
	        if (filetype==2) irplot([channel1,channel2],sampling_rate,title) ;return; end;
	        if (filetype==3)  logplot([channel1,channel2],sampling_rate,title,'gr',min_frq,max_frq); end;
    	   end; 

 	   if (channels_per_dir==1) 
		channel=mtlrch(catalogue,index);
		if (filetype==2) irplot(channel,sampling_rate,title); end;
		if (filetype==3) logplot(channel,sampling_rate,title,'',min_frq,max_frq); end;
   	   end;

	   axis([min_frq max_frq min_abs-max_abs 5])
           str_ele=int2str_(fix(elevation));
           str_azi=int2str_(fix(azimuth));

           graph_name= sprintf('lm%s%s.ps',str_ele,str_azi)
           eval(['print -dps ' graph_name]);
           close(clf);
           
        end;
    

end;

%if (max_frq > stop) max_frq=stop ;end;


if (nargin == 3)

	direction_matrix=mtlrdir(catalogue);
	dist=sphdist([elevation,azimuth],direction_matrix);
	if any(dist==0) 
   		index=find(dist==0);
	else 
  		[min_dist,index]= min(dist);
	end;
elevation=direction_matrix(1,index);
azimuth=direction_matrix(2,index);

title=sprintf('Catalogue : %s, Elevation: %i , Azimuth: %i', catalogue, elevation,azimuth);

end;

if (nargin == 2)  
	index=elevation; 
	title=sprintf('Catalogue : %s, Channel: %i', catalogue, index);
end;

if (index > n_channels) error('MTLPLOT: index > n_channel'); return;end;


if ((channels_per_dir==1) | (nargin ==2))
	channel=mtlrch(catalogue,index);
	if (filetype==2) irplot(channel,sampling_rate,title); grcntrl;return; end;
	if (filetype==3) 
	   	magplot(channel,sampling_rate,title,'',min_frq,max_frq);grcntrl; 
  		phasplot(channel,sampling_rate,title,'',min_frq,max_frq);grcntrl;
   		return; 
        end;
end;

if (channels_per_dir==2)
	channel1=mtlrch(catalogue,index*2-1);
	channel2=mtlrch(catalogue,index*2);
	
	if (filetype==2) irplot([channel1,channel2],sampling_rate,title); grcntrl;return; end;

	if (filetype==3) 
	   magplot([channel1,channel2],sampling_rate,title,'',min_frq,max_frq);grcntrl; 
   	   phasplot([channel1,channel2],sampling_rate,title,'',min_frq,max_frq);grcntrl;
   	   return; 
	end;

end;

return;



