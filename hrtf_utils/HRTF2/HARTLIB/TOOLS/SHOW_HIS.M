% SHOW_HIST Matlab script fuer die Darstellung  der Packet_Histogramm Daten
%
% (c) Christian Jung

function  y = show_hist(fname)

 eval(['load ',fname, '.dat' ])
 output = eval(fname);
 figure;
 groesse = size(output,1);
 groessestr = int2str(groesse);
 hist(output(:,3),400)
 title([fname, '.dat ,  Gesamtzahl der Versuche : ', groessestr ])
 xlabel('byte / s')
 ylabel('Anzahl')
 grcntrl;

