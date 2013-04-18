function out = mirror(n_lines,in,kind)

out = [];
if nargin==2
  out(1:n_lines) = in;
  for i=1:n_lines-2
   out(n_lines+i) = in(n_lines-i);
  end;
elseif nargin==3
  out(1:n_lines) = in;
  for i=1:n_lines-2
   out(n_lines+i) = conj( in(n_lines-i) );
  end;
end;