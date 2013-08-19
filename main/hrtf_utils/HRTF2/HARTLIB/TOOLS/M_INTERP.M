function yi = m_interp(x,y,xi,valid)

%function yi = m_interp(x,y,xi, [valid])
%
%x 	: original x-axis
%y 	: original y-data
%xi	: matrix with new x-axis
%valid	: array indicating number of valid rows per column in xi

% Klaus Hartung
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 7.2.96

if (nargin < 3) help m_interp; error('Some argments are missing'); end;

[m_max,n_max]= size(xi);

if (nargin < 4) valid = ones(n_max)*m_max; end;

yi = zeros(m_max,n_max);

for n=1:n_max
	for m=1:valid(n)
		yi(m,n)=interp1(x,y,xi(m,n));
	end
end 

return
