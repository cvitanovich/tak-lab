function frq_axis=mtlrfrq(filename)

% frq_axis=mtlrfrq(filename)
% 
% MTLrfrq reads frequency axis from MTLSIG catalogue (non-linear frequency axis)
%
%         frq_axis: array with frequency lines
%         filename: string with filename
% 

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 24.4.1995




[filetype,info_blocks,...
 n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,...
 comment1, comment2] = mtlrh(filename);

if ( (filetype >= 7) & (filetype <= 11) )

	[fid,message] = fopen(filename,'r','ieee-le');
	if (fid==-1)          % file can't be opened
  		disp(message);
  		frq_axis=-1;
  		return;
	end;


	status=fseek(fid,512+n_directions*4*2,'bof');
	if (status<0) error('ERROR: mtlrfrq,  fseek to frq. axis not succesful'); delay=-1; return; end;

	frq_axis=fread(fid,n_lines,'float32');

	fclose(fid);
	if (status<0) error('ERROR: mtlfrq,  fclose not succesful' ); frq_axis=-1; return; end;

	return;

else
error('ERROR: mtlfrq,  wrong filetype' ); 	
end;
