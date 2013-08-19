% Globals_GUI
% initializes all global structures for system II experiments
% using GUI

global GUI
GUI = struct(...
   'paradigm', 1,...
   'recordmode', 0,...
   'cursormode', 2,...
   'spaceres', 'max',...
   'hrtf_locations', [],...
   'locations1', [],...
   'locations2', [],...
   'mode', [],...
   'UseLastLocations_flag', 0);   
