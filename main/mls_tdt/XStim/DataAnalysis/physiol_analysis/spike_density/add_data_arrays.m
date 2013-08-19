clear
load d:\mlspezio\matlab\save\Neuron_24
%Test numbers
test_numbers = [
%[Thr  SR ABI ITD ABI/F TIF BIF IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12 TS1 TS2 TS3 TS4 TS5 TS6 TS7 TS8 TS9 TS10 TS11 TS12 ]
  -57  -1  -1  65    -1  66  68  70  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  71  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -62  -1  74  75    -1  76  77  79  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  80  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1  -1  85    -1  87  88  89  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  90  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -72  -1  -1  91    -1  92  93  94  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  95  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -80  -1  96  97    -1  98 100 101  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -72  -1 103 104    -1 105 106 107  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 108 109    -1 110 111 112  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 113  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 121 122    -1 123 124 125  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 126  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1  -1 127    -1 129 130 131  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 132  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -73  -1 133 134    -1 135 136 137  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 138  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -73  -1 141 142    -1 143 144 145  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 146  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 152 153    -1 162 158 159  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 161  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 165 166    -1 167 168 169  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 170  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -75  -1 171 172    -1 173 174 175  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 176  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 206 208    -1 209 210 211  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 207 212 213  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70 220 214 216   223 217 218 219  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 221  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 228 230   234 232 233 231  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 229  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 242 248    -1 249 250 247  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1 246  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 264 265    -1 274  -1 267 268  -1 269  -1 270 271 272 273   -1   -1   -1 266  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 289 290    -1 300  -1 293 294  -1 295  -1 296 297 298 299   -1   -1   -1 292  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -60  -1 318 319    -1  -1  -1 321 322  -1 323  -1  -1  -1  -1  -1   -1   -1   -1 320  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 326 327    -1 332  -1 329 330  -1 331 333 334 335  -1  -1   -1   -1   -1 328  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 337 338    -1 343  -1 340 341  -1 342  -1  -1  -1  -1  -1   -1   -1   -1 339  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 347 348    -1 356  -1 351 352  -1 353  -1  -1 355 358 357   -1   -1   -1 349  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
];


bins = [0:300];
cache_files_directory = 'd:\mlspezio\matlab\scripts\physiol_analysis\899analysis\899_cache';
for cell_num = 1:size(test_numbers,1)
%for cell_num = 1:1
%Tonal
   if(any(strcmp(fieldnames(Neuron{cell_num}),'tif_meansurf')));
   data_file_wildcard = sprintf('899rSITE*T%03g.MAT', test_numbers(cell_num,6));
   cachedir = eval(['dir(''' cache_files_directory '\' data_file_wildcard ''')']);
   if (size(cachedir,1)>1) error('Bizarre Error: More than one cache file matches this test number.'); end;
   if (size(cachedir,1)~=0) 
      data_file = cachedir(1).name;
      data_file = [cache_files_directory '\' data_file];
      disp(['Loading cached data from file: ' data_file])
		eval(['load ' data_file ';']);
   else
      mean_resp_surf = -1;
      disp('Error: Cached file doesnt exist, run proc_raw first.');
      return;
   end;
   
   Neuron{cell_num}.tif_dataarray = data_array;
	end
   
   %BP
   if(any(strcmp(fieldnames(Neuron{cell_num}),'bif_meansurf')));
   clear data_array
   data_file_wildcard = sprintf('899rSITE*T%03g.MAT', test_numbers(cell_num,7));
   cachedir = eval(['dir(''' cache_files_directory '\' data_file_wildcard ''')']);
   if (size(cachedir,1)>1) error('Bizarre Error: More than one cache file matches this test number.'); end;
   if (size(cachedir,1)~=0) 
      data_file = cachedir(1).name;
      data_file = [cache_files_directory '\' data_file];
      disp(['Loading cached data from file: ' data_file])
      eval(['load ' data_file ';']);
   else
      mean_resp_surf = -1;
      disp('Error: Cached file doesnt exist, run proc_raw first.');
      return;
   end;
   
	Neuron{cell_num}.bif_dataarray = data_array;
   end
   
   %ILDAlone
   if(any(strcmp(fieldnames(Neuron{cell_num}),'ia_meansurf')));
   clear data_array
   data_file_wildcard = sprintf('899rSITE*T%03g.MAT', test_numbers(cell_num,8));
   cachedir = eval(['dir(''' cache_files_directory '\' data_file_wildcard ''')']);
   if (size(cachedir,1)>1) error('Bizarre Error: More than one cache file matches this test number.'); end;
   if (size(cachedir,1)~=0) 
      data_file = cachedir(1).name;
      data_file = [cache_files_directory '\' data_file];
      disp(['Loading cached data from file: ' data_file])
      eval(['load ' data_file ';']);
   else
      mean_resp_surf = -1;
      disp('Error: Cached file doesnt exist, run proc_raw first.');
      return;
   end;
   
   Neuron{cell_num}.ia_dataarray = data_array;
	end
   
end
   
   
