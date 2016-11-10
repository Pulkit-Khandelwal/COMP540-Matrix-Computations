%calculating the given rotation parameters c and s where c computes the cos
% and s computes the sine of the theta. Theta is determined as below

function [c,s] = givensrotation(a,b)
  if b == 0
    c = 1;
    s = 0;
  else
      r = a / b;
      s = 1 / sqrt(1 + r^2);
      c = s*r;
    
  end

end