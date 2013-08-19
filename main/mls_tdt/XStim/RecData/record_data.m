function [] = record_data(params,datamatrix)
%Spike Data Recording Routine

datafile = [params.datadir params.datafile];
load(datafile);
if(isstruct(DATA))
   temp.Site = DATA.Site;
   clear DATA
   DATA = temp;
end
DATA.Site{params.recording_site}.test{params.testnum}.params = params;
DATA.Site{params.recording_site}.test{params.testnum}.datamatrix = datamatrix;
save(datafile,'DATA');
clear DATA temp

return;