%Add a neuron to the Neuron cell array
%This is meant for data collected the *new* way, not for data collected
%under David's DOS system

clear; close all

%Order of test_numbers:
%901L2 - Sites 7,8,9,11
%913R1 - Sites 2,3,8,9
%913R3 - Sites 3
%894L1 - Sites 3,7

%Bird Number:
%bird_number = 901
bird_number = 913;
%bird_number = 894;

%Side of Brain:
%side_of_brain = 'L';
side_of_brain = 'R';

%Session Number:
session_num = 1;

%Data Sites:
%data_sites = [7 8 9 11]; %901L2
data_sites = [2 3 8 9]; %913R1
%data_sites = [3]; %913R3
%data_sites = [3 7]; %894L1

%Threshholds:
%threshholds = [-60 -72 -80 -75]; %901L2
threshholds = [-72 -73 -80 -77]; %913R1
%threshholds = [-60]; %913R3
%threshholds = [-72 -80]; %894L1

test_numbers = [
%[SR ABL ITD ABI/F TIF TSIF GIF GSIF BBIF FD1 FD2 FD3 FD4 IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12 TS1 TS2 TS3 TS4 TS5 TS6 TS7 TS8 TS9 TS10 TS11 TS12 ]
%  -1   2   3    -1   8    7   6   -1   -1  -1  -1  -1  -1   5  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   4  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
%  -1   2   4    -1  -1    8   7   -1   -1  -1  -1  -1  -1   6  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   5  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
%  -1   2   3    -1   8    7   6    8   -1  -1  -1  -1  -1   5  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   4  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
%  -1   2   3    -1  -1    6   5    7   -1  -1  -1  -1  -1   4  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   8  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -1   3   5    -1  -1    6   8   10   11  -1  -1  -1  -1   6  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   2  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -1   2   5    -1  -1    6   7    9    8  11  12  -1  -1   6  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   4  14  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -1   2   4    -1  -1    6   7   13   10  11  14  -1  -1   6  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   5  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -1   3   5    -1  -1    6   8   -1   10  11  12  13  -1   7  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   6  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
 % -1   8   9    -1  -1   -1  12   -1   14  15  16  17  -1  11  18  19  -1  -1  -1  -1  -1  -1   -1   -1   -1  10  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
 % -1   3   4    -1  -1   -1   6   -1    7   8   9  -1  -1   5  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1   2  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
 % -1   2   3    -1  -1   -1   6   -1    9  10  11  12  13   5  15  16  17  18  -1  -1  -1  -1   -1   -1   -1   4  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
];


%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
data_file = ['DATA' num2str(bird_number) side_of_brain num2str(session_num)];
%hrtf_file       = 'd:\mlspezio\owl_hrtfdata\899\out9be';
%hrtf_file       = 'd:\mlspezio\owl_hrtfdata\901\901flat_freqdomAD';
hrtf_file       = 'd:\mlspezio\owl_hrtfdata\913\913flat_freqdomAD';
%hrtf_file       = 'd:\mlspezio\owl_hrtfdata\913\913f_fdAD';
%hrtf_file       = 'd:\mlspezio\owl_hrtfdata\894\894f_fdAD';
hrtf_file2		  = '';
get_HRTF = 0;
get_ITDmatrix = 0;
HRTFAll = [];

if(exist('d:\mlspezio\matlab\save\HRTFAll_2.mat'))
   eval(['load d:\mlspezio\matlab\save\HRTFAll_2'])
end
%Get HRTF info if desired
if(get_HRTF == 1)
   [ILD_matrix,ABL_matrix,ITD_matrix,HRTFinfo] = ...
      get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);
   numofHRTF = length(HRTFAll) + 1;
   HRTFAll{numofHRTF}.bird_number = bird_number;
   HRTFAll{numofHRTF}.HRTFinfo = HRTFinfo;
   HRTFAll{numofHRTF}.ILD_matrix = ILD_matrix;
   HRTFAll{numofHRTF}.ABL_matrix = ABL_matrix;
end
eval(['save d:\mlspezio\matlab\save\HRTFAll_' num2str(length(HRTFAll)) ' HRTFAll']);
ILD_matrix = HRTFAll{2}.ILD_matrix;
HRTFinfo = HRTFAll{2}.HRTFinfo;

%Begin to add the neuron structure
num_sites_to_add = length(data_sites);

%Load old Neuron structure
load 'd:\mlspezio\matlab\save\Neuron_28'

%Load data to be added
eval(['load ' data_directory data_file]);

incr = 0;
plotflag = 1;
start = length(Neuron) + 1;
for cell_num = start:start+num_sites_to_add - 1;
   incr = incr + 1;
   Neuron{cell_num} = ...
      make_neuronstruct2(DATA.Site{data_sites(incr)},...
      bird_number,...
      ILD_matrix,...
      HRTFinfo,...
      threshholds(incr),...
      test_numbers(incr,:),...
      plotflag);
end

eval(['save d:\mlspezio\matlab\save\Neuron_' num2str(length(Neuron)) ' Neuron']);
