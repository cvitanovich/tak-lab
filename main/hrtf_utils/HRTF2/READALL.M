% READALL  script to read a set of hrtf's from one recording session

session = '802D_'
data_path = 'hrtfdata/'
file_suffix = '.SC2'

start_index = 1;
end_index = 616;
%spectrum1 = []; spectrum2 = [];
%azim_all = []; elev_all = [];

for i = start_index:end_index

   % get number as string and pad with zeros
   num_string = int2str(i);
   while size(num_string)<3
      num_string = ['0',num_string];
   end;

   filename = [data_path, session, num_string, file_suffix]
   [spec1, spec2, az, elev, sr, sf, ef] = readhrtf(filename);

   % protect against misread
   if spec1==[]
      spec1 = zeros(1,1024);
   end;
   if spec2==[]
      spec2 = zeros(1,1024);
   end;
   if az == []; az = NaN; end;
   if elev == []; elev = NaN; end;

   spectrum1 = [spectrum1; spec1];
   spectrum2 = [spectrum2; spec2];
   azim_all = [azim_all; az]; elev_all = [elev_all elev];
   
end;
