function mtlplot(catalogue, elevation, azimuth)

% mtlplot(catalogue, elevation, azimuth)
%
% mtlplot(catalogue, channel)
%
% MTLPLOT plots the impulse response or transfer function of the selected
% direction or channel

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995

min_frq=100;
max_frq=18000;

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
start,stop,n_directions,comment1, comment2] = mtlrh(catalogue);

if (max_frq > stop) max_frq=stop ;end;

binaural=0;
if ((n_channels/n_directions)==2) binaural=1;

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


if ((binaural==0) | (nargin ==2))
      channel=mtlrch(catalogue,index);
      if (file_type==2) irplot(channel,sampling_rate,title); grcntrl;return; end;
      if (file_type==3) 
           magplot(channel,sampling_rate,title,'',min_frq,max_frq);grcntrl; 
           phasplot(channel,sampling_rate,title,'',min_frq,max_frq);grcntrl;
           return; 
        end;
end;

if (binaural==1)
      channel1=mtlrch(catalogue,index*2-1);
      channel2=mtlrch(catalogue,index*2);
      
      if (file_type==2) irplot([channel1,channel2],sampling_rate,title); grcntrl;return; end;

      if (file_type==3) 
         magplot([channel1,channel2],sampling_rate,title,'',min_frq,max_frq);grcntrl; 
         phasplot([channel1,channel2],sampling_rate,title,'',min_frq,max_frq);grcntrl;
         return; 
      end;

end;

return;
end;