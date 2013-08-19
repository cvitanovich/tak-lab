

dir = [35   -25
    30   -20
    25   -15
    20   -10
    15    -5
    10     0
    5     5
    0    10
    -5    15
    -10    20
    -15    25
    -20    30
    -25    35
    30   -30
    25   -25
    20   -20
    15   -15
    10   -10
    5    -5
    0     0
    -5     5
    -10    10
    -15    15
    -20    20
    -25    25
    -30    30
    25   -35
    20   -30
    15   -25
    10   -20
    5   -15
    0   -10
    -5    -5
    -10     0
    -15     5
    -20    10
    -25    15
    -30    20
    -35    25
    20   -40
    15   -35
    10   -30
    5   -25
    0   -20
    -5   -15
    -10   -10
    -15    -5
    -20     0
    -25     5
    -30    10
    -35    15
    -40    20
    15   -45
    10   -40
    5   -35
    0   -30
    -5   -25
    -10   -20
    -15   -15
    -20   -10
    -25    -5
    -30     0
    -35     5
    -40    10
    -45    15
    10   -50
    5   -45
    0   -40
    -5   -35
    -10   -30
    -15   -25
    -20   -20
    -25   -15
    -30   -10
    -35    -5
    -40     0
    -45     5
    -50    10
    5   -55
    0   -50
    -5   -45
    -10   -40
    -15   -35
    -20   -30
    -25   -25
    -30   -20
    -35   -15
    -40   -10
    -45    -5
    -50     0
    -55     5]';

dir(1,:) = dir(1,:)+15;
dir(2,:) = dir(2,:)-15;

[n nLocs] = size(dir);

for iLoc = 1:nLocs
    stim_dir = 'c:\kip_overflow\2sources\cues\ABL_90Hz\';
    el = dir(1,iLoc);   az = dir(2,iLoc);
    [Lnoi1, Rnoi1] = read_mls_stim(el,az, stim_dir);
    stim_dir = 'c:\kip_overflow\2sources\cues\ITA\';
    [Lnoi2, Rnoi2] = read_mls_stim(el,az, stim_dir);
    Lnoi2 = Lnoi2(128:36762-127);
    Rnoi2 = Rnoi2(128:36762-127);

    Lnoi = (Lnoi1 + Lnoi2)/2;
    Rnoi = (Rnoi1 + Rnoi2)/2;
    stim_dir = 'c:\kip_overflow\2sources\cues\Combined_55_75_90\';
    
    fname = [stim_dir 'stimbuf_left_' num2str(el) '_' num2str(az)];
    eval(['fid = fopen(''' fname ''', ''a'');'])
    fwrite(fid, Lnoi, 'short');
    
    fname = [stim_dir 'stimbuf_right_' num2str(el) '_' num2str(az)];
    eval(['fid = fopen(''' fname ''', ''a'');'])
    fwrite(fid, Rnoi, 'short');
    fclose(fid);
    clear L* R* stim_dir
end