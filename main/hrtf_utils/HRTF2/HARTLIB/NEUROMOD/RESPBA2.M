function resp_ba2(listname)

% RESP_BA2(listname)
%respfname='ca060101'
%

old_animal_number='';

list=readlist(listname);
[n_list,x] = size(list);
disp(['Number of recordings: ' int2str(n_list)])
for n=1:n_list
	respfname=list(n,:)
	T=readsmr(['/HRTF/cavia/smr/' respfname]);
	animal_number=respfname(1:4)
	penetration_number=respfname(5:6)
	recording_number=respfname(7:8)
	subst=[animal_number(1:2) int2str(str2num(animal_number(3:4)))];
		
	direc=['/HRTF/cavia/' subst]; 
	if (strcmp(old_animal_number,animal_number)==0)
		action=['load /HRTF/cavia/' subst '/' subst  'eq']
		
		disp('Loading  ');
		eval(action);
		action1=['load /HRTF/cavia/' subst '/' subst  'c']
		eval(action1);
	
		HRIR=HRIR_LS;
		clear HRTF;
		clear HRIR_LS;
		HRTF=rfft(HRIR,1024);
		clear HRIR;
		str=['/HRTF/cavia/' subst '/ITDILD.mat']
		fid=fopen(str,'r')
		if (fid >0)
			fclose(fid)
			eval(['load /HRTF/cavia/' subst '/ITDILD.mat']) 		
		 else	
			ana
%			save
			action=['save /HRTF/cavia/' subst '/ITDILD  ILD ITD LEVEL b_array']
			eval(action);
		 end;
		 old_animal_number=animal_number;
	end;

	recordname=strtok(respfname,'.');
	close all;
	anawdh(recordname,ILD,ITD,LEVEL,T,MAP,b_array,direc);
end
