function resp_ba(listname)

% RESP_BA batch for calculation of interaural and monaural parameter and weighting factors
%
% RESP_BA(list)
%
% list string with name of list

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% May 1996

%respfname='ca060101'
%

old_animal_number='';

list=readlist(listname);
[n_list,x] = size(list);
disp(['Number of recordings: ' int2str(n_list)])
for n=1:n_list
	respfname=lower(list(n,:))
	T=readsmr(['/home/aea/hartung/hrtf/cavia/smr/' respfname]);
	animal_number=respfname(1:4)
	penetration_number=respfname(5:6)
	recording_number=respfname(7:8)
	subst=[animal_number(1:2) int2str(str2num(animal_number(3:4)))];
		
	direc=['/home/aea/hartung/hrtf/cavia/' subst]; 
	if (strcmp(old_animal_number,animal_number)==0)
		action=['load ' direc '/' subst  'eq']
		
		disp('Loading  ');
		eval(action);
		action1=['load ' direc '/' subst  'c']
		eval(action1);
	
		HRIR=HRIR_LS;
		clear HRTF;
		clear HRIR_LS;
		HRTF=rfft(HRIR,1024);
		clear HRIR;
		str=[direc '/ITDILD.mat']
		fid=fopen(str,'r')
		if (fid >0)
			fclose(fid)
			eval(['load ' direc '/ITDILD.mat']) 		
		 else	
			ana
%			save
			action=['save ' direc '/ITDILD  ILD ITD LEVEL b_array']
			eval(action);
		 end;
		 old_animal_number=animal_number
	end;

	recordname=strtok(respfname,'.');
	close all;
	anaresp(recordname,ILD,ITD,LEVEL,T,MAP,b_array,direc);
end
