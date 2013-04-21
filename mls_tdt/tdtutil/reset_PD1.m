function [status] = reset_PD1(status);

% checks TDT.PD1status and then resets S232 and PD1, PA4

global TDT

if TDT.S232_status
    S232('PA4mute',1);
    S232('PA4mute',2);
    S232('PD1stop',1);
    S232('PD1clear',1);
    S232('trash');
    S232('dropall');
    S232('S2close') % close application and release AP2 and XBUS locks
    TDT.S232_status = 0;
end
status = 0;