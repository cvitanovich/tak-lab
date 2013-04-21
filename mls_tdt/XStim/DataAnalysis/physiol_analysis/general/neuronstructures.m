%Script to use the function make_neuronstructure to take data
%and build a structure called a "Neuron"
%Look at the function header for a list of Neuron properties
%Data needs to have been gone through the raw processing already
%Test_Numbers format:
%[Thr SR ABI ITD ABI/F TIF BIF IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12 TS1 TS2 TS3 TS4 TS5 TS6 TS7 TS8 TS9 TS10 TS11 TS12 ]

clear; close all

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
];

reps = zeros(8,size(test_numbers,2),size(test_numbers,1));
tests = [98 112 123 129 169 170 211 212 213];
for x = 1:length(tests)
   switch tests(x)
   case 98,
      [row,col] = find(test_numbers == 98);
      reps(:,col,row) = 1:8;
   case 112,
      [row,col] = find(test_numbers == 112);
      reps(:,col,row) = 1:8;
   case 123,
      [row,col] = find(test_numbers == 123);
      reps(:,col,row) = [1:6 zeros(1,2)];
   case 129,
      [row,col] = find(test_numbers == 129);
      reps(:,col,row) = 1:8;
   case 169,
      [row,col] = find(test_numbers == 169);
      reps(:,col,row) = [1:5 zeros(1,3)];
   case 170,
      [row,col] = find(test_numbers == 170);
      reps(:,col,row) = [1:3 zeros(1,5)];
   case 211,
      [row,col] = find(test_numbers == 211);
      reps(:,col,row) = [1:5 zeros(1,3)];
   case 212,
      [row,col] = find(test_numbers == 212);
      reps(:,col,row) = [1:3 zeros(1,5)];
   case 213,
      [row,col] = find(test_numbers == 213);
      reps(:,col,row) = [1:3 zeros(1,5)];      
   end
end

%Parameters
bird_number = 899;
side_of_brain = 'r';

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file       = ['out9be'];
hrtf_file2		 = ['out29be'];
get_HRTF = 0;
get_ITDmatrix = 0;

[ILD_matrix,ABI_matrix,ITD_matrix,HRTFinfo] = get_HRTFinfo(bird_number,hrtf_file,hrtf_file2,get_HRTF,get_ITDmatrix);

if(0)
for cell_num = 1:size(test_numbers,1)
   for test_number = 1:size(test_numbers,2)
      if(test_numbers(cell_num,test_number) ~= -1)
         disp(['Cell ' num2str(cell_num) ', Test ' num2str(test_numbers(cell_num,test_number)) ':'])
         pause(1)
         [meansurf, stdsurf, dim2vals, dim1vals, testpars, spont_spikes, spont_dur, nincl_reps] = ...
            proc_test899(bird_number, side_of_brain, test_numbers(cell_num,test_number), 1, 0);
         pause(2)
      end
   end
end
end

for cell_num = 1:size(test_numbers,1)
   Neuron{cell_num} = ...
      make_neuronstruct(bird_number,side_of_brain,...
      test_numbers(cell_num,:),...
      ILD_matrix,HRTFinfo,reps(:,:,cell_num));
end

save 'd:\mlspezio\matlab\save\Neuron_base' Neuron