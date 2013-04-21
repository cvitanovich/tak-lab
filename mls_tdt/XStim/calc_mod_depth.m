function [depth] = calc_mod_depth(loud,soft)

%function [depth] = calc_mod_depth(loud,soft)
% input two dB levels
% output percent depth of modulatiom

if soft > loud
    error('ERROR    function [depth] = calc_mod_depth(loud,soft)');
end
depth = 100 - 100/ 10^((loud-soft)/20);
disp(['Percent Modulation Depth : ' num2str(depth)]);