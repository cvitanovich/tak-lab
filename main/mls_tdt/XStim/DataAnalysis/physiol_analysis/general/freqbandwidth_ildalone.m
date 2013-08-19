%Script to analyze the impact on ILDAlone Response Surface (IRS)
%of widening frequency bandwidth of a sound around its frequency of "most sensitivity"
clear;close all

%Tests for analysis
%Each row is a different cell
%Numbers in the following order: [1 kHz, 3kHz, 7kHz, 9kHz]
test_numbers = [
   219 224 225 226 -1
   231 235 236 237 238];
test_types = str2mat('0-13 kHz','1 kHz','3 kHz', '7 kHz', '9 kHz');
test_bands = [
   1 3 7 0
   1 3 7 9];

%Set parameters
bird_number = 899;
side_of_brain = 'r';
get_hrtf = 0;
fms = [5800 6275];

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];

%Get information from bird's HRTF catalogue - for use in predicting ILD Alone Surfaces
%Header information
[filetype,...
info_blocks,...
n_channels,n_lines, sampling_rate,...
first_line,last_line,num_locations,...
comment1, comment2] = mtlrh([hrtf_directory hrtf_file]);
%Location information
temp = mtlrdir([hrtf_directory hrtf_file]);
location_matrix = sph2dbl(temp);
itddir = location_matrix;
clear temp

hrtf_freqs = [first_line:(last_line-first_line)/(n_lines-1):last_line]';

for cell_num = 1:size(test_numbers,1)
   for test_num = 1:length(test_numbers(cell_num,:))
      if(test_numbers(cell_num,test_num) ~= -1)
   
%1. Get the mean response surface of the ILD Alone data
[ILDAlone_meansurf, ILDAlone_stdsurf, ILDAlone_dim2vals, ILDAlone_dim1vals, testpars, spont_spikes, spont_dur, nincl_reps, locs] = ...
   proc_test899(bird_number, side_of_brain, test_numbers(cell_num,test_num), 1, 0);
ILDAlone_meanmat(test_num,:) = ILDAlone_meansurf;
ILDAlone_stddmat(test_num,:) = ILDAlone_stdsurf;
[ILDAlone_azi,ILDAlone_ele,temp] = array2diamond(ILDAlone_meansurf,locs');
%Interpolate the missing values in the ILDAlone measurement
	[AZ EL] = meshgrid(ILDAlone_azi, ILDAlone_ele);

	% generate mask for missing points
	missmask = NaN*ones(size(temp));
	
	i = 1;
	for az = -90:5:90;
	  for el = -90+abs(az):5:90-abs(az)
	    if (~(locs(:,1)==el & locs(:,2)==az))
	       missmask(AZ==az & EL==el) = 1;
	    else
	       missmask(AZ==az & EL==el) = 0;
	    end; 
	  end;
	end;
	
	% replace all valid locations presently containing NaN's with zeros
	ri = (missmask==1);
	temp(ri) = zeros(size(temp(ri)));
	% generate interpolation function - uses surrounding 4 squares to generate missing square
   intfun = [[0 1 0]; [1 0 1]; [0 1 0]];
	temp(find(isnan(temp))) = zeros(size(find(isnan(temp))));
	intval_old = conv2(temp, intfun,'same')/4;
	intval = intval_old.*missmask;
	
   ILDAlone_diamond = intval + temp;
   ILDAlone_diammat(:,:,test_num) = ILDAlone_diamond;
   ILDAlone_diammat(:,:,test_num) = ILDAlone_diammat(:,:,test_num)./max(max(ILDAlone_diammat(:,:,1)));
   
end %end if statement
end %end loop over test numbers

figure
for test_num = 1:size(ILDAlone_meanmat,1)
   subplot(3,2,test_num)
   plotdiam(ILDAlone_azi,ILDAlone_ele,ILDAlone_diammat(:,:,test_num));
   set(gca,'Color','none');
   if(test_num == 1)
      title(['IRS for BW of ' deblank(test_types(test_num,:))])
   end
   if(test_num > 1)
      temp = corrcoef(ILDAlone_meanmat(1,:),ILDAlone_meanmat(test_num,:));
      cc(cell_num,test_num) = temp(1,2);
      title(['IRS for BW of ' deblank(test_types(test_num,:)) ' around '...
            num2str(fms(cell_num)) ', r = ' num2str(cc(cell_num,test_num))])
   end
   colorbar
end

clear ILDAlone_meanmat ILDAlone_stddmat ILDAlone_diammat

end %end loop over cells

figure
hold on
colors = ['g','m'];
symbols = ['*','s'];
for cell_num = 1:size(test_numbers,1)
   plot(dezero(test_bands(cell_num,:)),dezero(cc(cell_num,2:5)),[colors(cell_num) symbols(cell_num) '-'])
   text(1,dezero(cc(cell_num,2)),...
      ['\leftarrow\fontname{times}FMS = ' num2str(fms(cell_num))],'FontSize',10,'Color','black')
end
xlabel('Bandwidth around FMS (kHz)')
ylabel('r between ILDAlone Response Surfaces')