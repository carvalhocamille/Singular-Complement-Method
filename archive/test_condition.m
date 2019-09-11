function [z] = test_condition(t)
f = @(t) t.^2;
if (t>1)
    z = f(t);
else
    z = -f(t);
    
end