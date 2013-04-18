function [m,p] = max_env(h)

%function [m,p] = max_env(h)
%
% MAX_ENV calculates the maximum value and its position of the envelop
% of an impusle response

env=etc(h);
[m,p]=max(env);

end;
