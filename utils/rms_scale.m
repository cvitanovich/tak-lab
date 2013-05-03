function x=rms_scale(x,target_rms)
% scales a vector to a desired rms
old_rms = norm(x)/sqrt(length(x));
x = x .* (target_rms/old_rms);