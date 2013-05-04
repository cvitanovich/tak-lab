function soundBufferSetup
global PDR

for i0=1:PDR.SOUNDS_num_carriers
    len_buf = PDR.buf_pts;
    len_stim = length(PDR.LEAD_sounds{i0});
    lead_stim=zeros(1,len_buf);
    lag_stim=lead_stim;
    pad = PDR.buf_pts-len_stim; %ceil(len_buf-len_stim);
    lead_stim = [zeros(1,2) PDR.LEAD_sounds{i0} zeros(1,pad)];
    lag_stim = [zeros(1,2) PDR.LAG_sounds{i0} zeros(1,pad)];
    PDR.LEAD_sounds{i0} = lead_stim(1:PDR.buf_pts);
    PDR.LAG_sounds{i0} = lag_stim(1:PDR.buf_pts);
end