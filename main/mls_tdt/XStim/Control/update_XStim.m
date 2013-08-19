% update_XStim control GUI
% should be called at the end of each setinfo and each exit

set(H.Xstim_silence_lead,'string',XStimParams.silence_lead);
set(H.Xstim_silence_trail,'string',XStimParams.silence_trail);
set(H.Xstim_curr_stimdur,'string', XStimParams.curr_stimdur);
set(H.Xstim_test_ISI,'string', XStimParams.test_ISI);
