% script to read files for Koppl's Fig 10
% koppl_10 
% each was measured using standard graph paper


fid = fopen('d:\kip\temp\koppl10a.txt');
x = fscanf(fid,'%f');
fclose all

j=1;
for i = 1:2:length(x)
   LOfreq(j) = x(i);
   LOdb(j) = x(i+1);
   j=j+1;
end
disp('finished reading and parsing koppl10a');

fid = fopen('d:\kip\temp\koppl10b.txt');
x = fscanf(fid,'%f');
fclose all

j=1;
for i = 1:2:length(x)
   HIfreq(j) = x(i);
   HIdb(j) = x(i+1);
   j=j+1;
end
disp('finished reading and parsing koppl10b');

fid = fopen('d:\kip\temp\koppl10c.txt');
x = fscanf(fid,'%f');
fclose all

j=1;
for i = 1:2:length(x)
   TOPfreq(j) = x(i);
   TOPdb(j) = x(i+1);
   j=j+1;
end
disp('finished reading and parsing koppl10b');

% convert to real numbers from graph-pad numbers
TOPdb = TOPdb/2.04;
HIdb = HIdb/1.57;
LOdb = LOdb/1.57;
TOPfreq = exp((2+2*TOPfreq/204)*log(10));   
HIfreq = exp((2+ 2*HIfreq/157)*log(10));     % now in log10 units
LOfreq = exp((2+ 2*LOfreq/157)*log(10));    % log10 units

disp('finished');

semilogx(HIfreq,HIdb,'r.')
hold on
semilogx(LOfreq,LOdb,'.')
semilogx(TOPfreq,TOPdb,'g.')
fclose all
