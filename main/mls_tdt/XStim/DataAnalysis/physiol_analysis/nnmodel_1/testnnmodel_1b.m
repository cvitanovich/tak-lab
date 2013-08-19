clear;close all
nntwarn off %needed to suppress warning messages regarding obsolete neural network functions

%Input needed for functions below
bird_number = 899;
side_of_brain = 'r';
test_numbers = [
   76 77 79 80 75 -62
   92 93 94 95 91 -72
   98 100 101 -1 97 -80
   105 106 107 -1 104 -72
   110 111 112 113 109 -75
   123 124 125 126 122 -75
   129 130 131 132 127 -65
   135 136 137 138 134 -73
   143 144 145 146 142 -73
   162 158 159 161 153 -75
   167 168 169 170 166 -70
   173 174 175 176 172 -75];
colormap_var = 'jet';
plotflag = 0;

%Directories & files
data_directory  = ['d:\mlspezio\physioldata\' num2str(bird_number) '\'];
cache_directory = ['d:\mlspezio\matlab\scripts\physiol_analysis\' num2str(bird_number) 'analysis\' ...
      num2str(bird_number) '_cache\'];
hrtf_directory  = ['d:\mlspezio\owl_hrtfdata\' num2str(bird_number) '\'];
hrtf_file = 'out9be';
get_hrtf = 0;

   figure(1)
   [net,New_BP_surface] = ...
      nnmodel_2(bird_number,side_of_brain,...
      test_numbers,hrtf_file,get_hrtf);

