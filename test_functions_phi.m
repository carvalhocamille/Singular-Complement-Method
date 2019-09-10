
[zeta, eta, etap, etapp] = cutoff(0.1, 0.6, 0.001);
phi  = pi/ 6;
Racine = 1.322404244;

Phi1  = @(t) sinh(Racine*(pi-t)) / sinh(Racine*(pi-phi/2)) ;
Phi2  = @(t) sinh(Racine*(t)) / sinh(Racine*(phi/2)) ;
Phi3  = @(t) -sinh(Racine*(pi+t)) / sinh(Racine*(pi-phi/2)) ;

Decay = @(r) r.^(1i * Racine) ;
DecayLog = @(z) exp(1i * Racine * z) ;


%% Create the singularity s(r,t) = r^lambda Phi(t)
Singu1 = @(r,t) Decay(r) .* Phi1(t);
Singu2 = @(r,t) Decay(r) .* Phi2(t);
Singu3 = @(r,t) Decay(r) .* Phi3(t);

% Display in the (log(r),t) coordinate system
f1 = @(z,t) real(DecayLog(z) .* Phi1(t));
f2 = @(z,t) real(DecayLog(z) .* Phi1(t));
f3 = @(z,t) real(DecayLog(z) .* Phi1(t));

figure;
fsurf(f1, [-10 0 phi/2 pi]);
xlabel('z');
ylabel('t');
figure;
fsurf(f2, [-10 0 -phi/2 phi/2]);
xlabel('z');
ylabel('t');
figure;
fsurf(f3, [-10 0 -pi -phi/2]);
xlabel('z');
ylabel('t');
%% Create other terms needed for computations
%Compute div(sigma Grad eta s)
Term = @(r) (1 + 2* 1i* Racine) * r .* etap(r) + r.^2 .* etapp - 2 * Racine^2 * eta(r);

DivSigma1 = @(r,t) sigma_d * (r.^(1i * Racine -2)) .* Phi1(t) .* Term(r);
DivSigma2 = @(r,t) sigma_m * (r.^(1i * Racine -2)) .* Phi2(t) .* Term(r);
DivSigma3 = @(r,t) sigma_d * (r.^(1i * Racine -2)) .* Phi3(t) .* Term(r);
