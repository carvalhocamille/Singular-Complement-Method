function [Beta_s,div_s_grad_s, s_grad_s1,Zsingu,sPhi,Phi] = singular_expression_radial(eta, etap, etapp, Racine,phi,sigma_m,sigma_d)
% Compute the radial componants o
%
%
%
%
div_s_grad_r  = @(r) (r*(1 + 2 * 1i * Racine) .* etap(r) + r.^2 .* etapp(r) - 2 * Racine^2 * eta(r)) .* r.^(1i * Racine -2);
grad_r = @(r) (r.* etap(r) + 1i * Racine * eta(r)) .* r.^(1i * Racine -1);
grad_t = @(r) eta(r) .* r.^(1i * Racine -1);


Phi1  = @(t) sinh(Racine*(pi-t)) / sinh(Racine*(pi-phi/2));
Phi2  = @(t) -sinh(Racine*(pi+t)) / sinh(Racine*(pi-phi/2));
Phi3  = @(t) sinh(Racine*(t)) / sinh(Racine*(phi/2));

Phip1  = @(t) -Racine*cosh(Racine*(pi-t)) / sinh(Racine*(pi-phi/2));
Phip2  = @(t) -Racine*cosh(Racine*(pi+t)) / sinh(Racine*(pi-phi/2));
Phip3  = @(t) Racine*cosh(Racine*(t)) / sinh(Racine*(phi/2));

[Phi] = @(t) (create_phi(Phi1, Phi2, Phi3, t, phi));
[sPhi] = @(t) (create_sphi(Phi1, Phi2, Phi3, t, phi, sigma_m, sigma_d));

[Phip] = @(t) (create_phip(Phip1, Phip2, Phip3, t, phi));
[sPhip] = @(t) (create_sphip(Phip1, Phip2, Phip3, t, phi, sigma_m, sigma_d));

%% Create expressions involving the singular part tobe computed via NCQ
div_s_grad_s = @(r,t) sPhi(t).* div_s_grad_r(r); % TEST sans sigma
s_grad_s1 = @(r,t) grad_r(r).* Phi(t) ;
s_grad_s2 = @(r,t) grad_t(r).* Phip(t) ;

Zsingu = @(r,t) eta(r).* r.^(1i*Racine) .* Phi(t);

%% Create exact terms involving the singular part only
Beta_s  = 0;

end

function [out] = create_phi(Phi1, Phi2, Phi3, t, phi)
if ( phi/2 <= t)
    out  = Phi1(t);
elseif  ( t <= -phi/2 )
    out   = Phi2(t);
else
    out  = Phi3(t);
end

end
    
function [out] = create_sphi(Phi1, Phi2, Phi3,t, phi, sigma_m, sigma_d)
if ( (phi/2 <= t) )
    out  = sigma_d * Phi1(t);
elseif ( t <= -phi/2  )
    out  = sigma_d * Phi2(t);
else
    out  = sigma_m *Phi3(t);
end

end
    
function [out] = create_phip(Phip1, Phip2, Phip3, t, phi)
if ( phi/2 <= t)
    out  = Phip1(t);
elseif  ( t<= -phi/2 )
    out   = Phip2(t);
else
    out  = Phip3(t);
end

end
    
function [out] = create_sphip(Phip1, Phip2, Phip3,t, phi, sigma_m, sigma_d)
if ( phi/2 <= t)
    out  = sigma_d * Phip1(t);
elseif ( t <= -phi/2)
    out  = sigma_d * Phip2(t);
else
    out  = sigma_m *Phip3(t);
end

end
    
    
    
    
    