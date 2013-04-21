function [out_sound, sound_specs] = amstim3(sample_rate, stimulus_length, amplitude_matrix, length_matrix)
% AMSTIM2(fn, sample rate, stimulus length, amplitude matrix, length matrix)
%
% Generate variance switching AM stimulus, where each segment is specifiable.
%
% fn: output file name.
%
% sample rate: the number of samples per second.
%
% stimulus length: the overall stimulus length in milliseconds.
%
% amplitude_matrix: the variance of square wave modulated noise switches between
%   				size(amplitude_vector, 2) segments in order, and after the last,
%					repeats with the first.  Segment N is constructed by multiplying
%					broadband noise by amplitude_vector(1, N).  If amplitude_vector(1, N)
%					is -1, then segment N is constructed by multiplying a random variable
%					drawn from a distrubution with param1 =  amplitude_vector(3, N) and
%					param2 = amplitude_vector(4, N).  The distribution is specified by the value
%					of amplitude_vector(2, N), where 1 = uniform and 2 = normal.  In the case
%					of uniform distributions, param1 = min and param2 = max, and in the case
%					of normal distributions, param1 = mean and param2 = standard deviation.
%
% length_matrix: same as amplitude_matrix, except for the lengths of the segments.
	
	count1 = 0;
	count3 = 1;
	while count1 <= (sample_rate * (stimulus_length / 1000))
		for count2 = 1:size(length_matrix, 2)
			
			if length_matrix(1, count2) ~= -1
				temp_length = length_matrix(1, count2);
			elseif length_matrix(2, count2) == 1
				temp_length = 0;
				while temp_length < 1
					temp_length = unifrnd(length_matrix(3, count2), length_matrix(4, count2));
				end
			else
				temp_length = 0;
				while temp_length < 1
					temp_length = normrnd(length_matrix(3, count2), length_matrix(4, count2));
				end
			end
				
			if amplitude_matrix(1, count2) ~= -1
				temp_amplitude = amplitude_matrix(1, count2);
			elseif amplitude_matrix(2, count2) == 1
				temp_amplitude = 0;
				while temp_amplitude <= 0
					temp_amplitude = unifrnd(amplitude_matrix(3, count2), amplitude_matrix(4, count2));
				end
			else
				temp_amplitude = 0;
				while temp_amplitude <= 0
					temp_amplitude = normrnd(amplitude_matrix(3, count2), amplitude_matrix(4, count2));
				end
			end
					
			temp_sound = MakeBBNoise(sample_rate, temp_length);
			
			temp_sound = temp_sound * temp_amplitude;
				
			if (count1 == 0) & (count2 == 1)
				out_sound = temp_sound;
			else
				out_sound = [out_sound temp_sound];
			end
			
            temp_amplitude
			sound_specs(count3, 1) = temp_amplitude;
			sound_specs(count3, 2) = temp_length;
			sound_specs(count3, 3) = length(out_sound) - length(temp_sound) + 1;
			sound_specs(count3, 4) = length(out_sound);
			count3 = count3 + 1;
		end
		
		count1 = length(out_sound);
	end
	
	if length(out_sound) > (sample_rate * (stimulus_length / 1000))
		out_sound(((sample_rate * (stimulus_length / 1000)) + 1):length(out_sound)) = [];
	end
end

