%Add a neuron to an existing Neuron cell array
%Script to use the function make_neuronstructure to take data
%and build a structure called a "Neuron"
%Look at the function header for a list of Neuron properties
%Data needs to have gone through the raw processing already
%Test_Numbers format:
%[Thr SR ABI ITD ABI/F TIF BIF IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12 TS1 TS2 TS3 TS4 TS5 TS6 TS7 TS8 TS9 TS10 TS11 TS12 ]

clear; close all

load 'd:\mlspezio\matlab\save\Neuron_20'
temp_cell = Neuron;
clear Neuron
Neuron = cell(1,length(temp_cell) + size(test_numbers,1));
for cell_num = 1:length(temp_cell)
   Neuron{cell_num} = temp_cell{cell_num};
end


test_numbers = [
%[Thr  SR ABI ITD ABI/F TIF BIF IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12 TS1 TS2 TS3 TS4 TS5 TS6 TS7 TS8 TS9 TS10 TS11 TS12 ]
% -70  -1 264 265    -1 274  -1 267 268  -1 269  -1 270 271 272 273   -1   -1   -1 266  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
% -65  -1 289 290    -1 300  -1 293 294  -1 295  -1 296 297 298 299   -1   -1   -1 292  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -60  -1 318 319    -1  -1  -1 321 322  -1 323  -1  -1  -1  -1  -1   -1   -1   -1 320  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -65  -1 326 327    -1 332  -1 329 330  -1 331 333 334 335  -1  -1   -1   -1   -1 328  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 337 338    -1 343  -1 340 341  -1 342  -1  -1  -1  -1  -1   -1   -1   -1 339  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
  -70  -1 347 348    -1 356  -1 351 352  -1 353  -1  -1 355 358 357   -1   -1   -1 349  -1  -1  -1  -1  -1  -1  -1  -1   -1   -1   -1
];


%Parameters
bird_number = 899;
side_of_brain = 'r';

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file       = [hrtf_directory 'out9be'];
hrtf_file2		 = [hrtf_directory 'out29be'];
get_HRTF = 1;
get_ITDmatrix = 0;

[ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);


reps = zeros(1,size(test_numbers,2),length(temp_cell)+size(test_numbers,1));

start = length(temp_cell) + 1;
for cell_num = start:(start + size(test_numbers,1) - 1)
   Neuron{cell_num} = ...
      make_neuronstruct(bird_number,side_of_brain,test_numbers(cell_num-start+1,:),ILD_matrix,HRTFinfo,reps(:,:,cell_num));
end

eval(['save d:\mlspezio\matlab\save\Neuron_' num2str(length(Neuron)) ' Neuron']);
