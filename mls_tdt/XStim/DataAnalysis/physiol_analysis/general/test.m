clear;close all

bird_number = 899;
side_of_brain = 'r';
test_numbers = [79];
hrtf_file = 'out9be';
get_hrtf = 0;

[optimal_ildf_surface] = optimal_ildf(bird_number,side_of_brain,test_numbers,hrtf_file,get_hrtf)