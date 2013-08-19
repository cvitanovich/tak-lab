%Exit_tdt: shuts down PD1 cleanly for good exit

S232('PD1stop',1);
S232('PD1clear',1);
S232('PA4mute',1);
S232('PA4mute',2);
S232('APunlock',0);
S232('XBunlock',0);
