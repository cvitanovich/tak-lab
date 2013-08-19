function status = exist2(S)

global GUI
global H
global XStimParams
global FN

warning('off');
X = [''' '''];
SS = [X(1) S X(3)];
pt = findstr(S,'.');
len = length(S);
if isempty(pt)
    pt = len +1;
end
S_name = S(1:pt-1);			
F_name = S(pt+1:len);
status = 0;
if isfield(eval(S_name), F_name)
	F = getfield(eval(S_name), F_name);
    if ~isempty(F);			
	    status = 1;
    end
end
warning('on');

