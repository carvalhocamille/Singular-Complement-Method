% cut-off function
delta = 0.1;
l = 0.6;
x = [delta: 0.001:l]; 
zeta = @(x) 1 - (1/2)*(1+erf((((x-delta)/l)*6-2.)/(sqrt(2)*0.3)));
zeta(delta)
zeta(l)
figure;
plot(x,zeta(x))