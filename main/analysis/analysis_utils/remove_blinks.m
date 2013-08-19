function A=remove_blinks(A)
%Plot fig to pick blink threshold levels
hblink = figure('Name','Remove A.trials with Blinks','NumberTitle','off');
set(gcf, 'Position', [0.1 0.1 0.4*A.scrn(3) 0.4*A.scrn(4)]);
plot (A.meandata, 'k'); grid; zoom on; set (gcf, 'pointer', 'crosshair');

% Added this to give a marker of sound onset (AC 6/27/2012):
A.startInt = A.bufs_before*(A.dec_pts-A.infopts)+round((A.dec_pts-A.infopts)*A.sndStart); % start integrating when sound starts in the buffer
line([A.startInt A.startInt],[min(min(A.meandata)) max(max(A.meandata))]);
axis tight
A.lolimit = input ('Lower cutoff =  ');
A.hilimit = input ('Upper cutoff =  ');
close(hblink);

%FILTER OUT ALL A.trials WITH BLINKS

for trial = 1:A.ntrials
    loblinks = min(find((A.meandata(:,trial)) <= A.lolimit));
    hiblinks = min(find((A.meandata(:,trial)) >= A.hilimit));
    loblinks = isempty (loblinks);
    hiblinks = isempty (hiblinks);
    loblink (trial) = loblinks;
    hiblink (trial) = hiblinks;
    
end

lonay = find (loblink ==1);
hinay = find (hiblink ==1);
A.nay = intersect(lonay, hinay);
A.nbltrials=(length(A.nay));

%cull out blink-free A.trials

A.mdata = zeros(A.bufs_total*(A.dec_pts-A.infopts),A.nbltrials)*NaN;

for n = 1:A.nbltrials;
    A.mdata(:,n) = A.meandata(:,A.nay(n));
end