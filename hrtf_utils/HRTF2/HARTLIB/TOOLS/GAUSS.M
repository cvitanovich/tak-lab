function result = te gauss(value)

difference = fix(value) - value;
if (difference~=0)
  result = fix(value) + 1;
else
  result = fix(value);
end;
