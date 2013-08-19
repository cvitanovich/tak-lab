function [filtered_buffer, circ_buffer] = circ_fir(circ_buffer,new_buffer,coefs)
    % takes inputs and updates a circular buffer with FIR filtered inputs
    % circ_buffer -- circular buffer 
    % (should be at least length of new_buffer + length of nTAPS)
    % new_buffer -- the new input buffer that gets filtered
    % coefs -- FIR coefficients
    % nTAPS -- number of coefficients
    % filtered_buffer -- the filtered output
    nTAPS=length(coefs);
    dur=length(new_buffer);
    tmp=NaN.*ones(1,(dur+nTAPS-1));
    circ_buffer=[circ_buffer((end-nTAPS+1):end) new_buffer]; % circular buffer for data (shift by duration of each data buffer)
    tmp=conv(coefs,circ_buffer);
    filtered_buffer=tmp(nTAPS+1:dur+nTAPS);