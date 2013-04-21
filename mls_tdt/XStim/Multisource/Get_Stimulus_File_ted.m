function [stim_from_file,stim,path] = Get_Stimulus_File_ted(prev_file, path)

[stim,path] = uigetfile([path '*.noi'],'Select Stimulus File');

if(path ~= 0)
   fid = fopen([path stim],'r');
   stim_from_file = fread(fid,inf,'float');
   fclose(fid);
else
   fid = fopen(prev_file,'r');
   stim_from_file = fread(fid,inf,'float');
   fclose(fid);
end

