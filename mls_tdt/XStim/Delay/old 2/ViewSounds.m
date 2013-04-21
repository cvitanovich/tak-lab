% ViewSounds 
% loads and displays TDT sounds

[FileName,PathName,FilterIndex] = uigetfile('*.*','Select TDT Sound');
if(FilterIndex==0), return, end

disp([PathName FileName]);

if(findstr('left', FileName) > 0)
    L_name = FileName;
    R_name = strrep(FileName,'left', 'right');
else
    L_name = strrep(FileName,'right', 'left');
    R_name = FileName;
end
%disp(L_name);
%disp(R_name);

fid=fopen([PathName L_name], 'r', 'l');
L_snd = fread(fid, inf, 'int16');

fid=fopen([PathName R_name], 'r', 'l');
R_snd = fread(fid, inf, 'int16');


if(findobj('Tag',FileName))
    close(FileName);
end
Stim.fig = figure('Units','pixels',...
  'Position',[50 500 500 300],...
  'Tag', FileName,...
  'Name',FileName,...
  'NumberTitle','off',...
  'Color',[1 1 1]);

%subplot(2,1,1);
Stim.L = plot(L_snd + 0.2*32767,...
    'Marker', 'none',...
    'LineStyle','-',...
    'Color', [1 0 0],...
    'LineWidth', 0.5);
%subplot(2,1,2);
hold on;
Stim.R = plot(R_snd - 0.2*32767,...
    'Marker', 'none',...
    'LineStyle','-',...
    'Color', [0 0 1],...
    'LineWidth', 0.5);
hold off;

