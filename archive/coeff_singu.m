function [Phi2theta] = coeff_singu(Numtri2,Coorneu,Refneu,Nbpt,contrast,RacineR)
% Formulas are the fomrulas (22), (23) in Bonnet-Ben Dhia, Carvalho, Chesnel, Ciarlet, JCP, 2016.
%% -------- Creation du vecteur Phi2theta ----------------------    
phi = pi/6;
RR = sqrt(Coorneu(:,1).*Coorneu(:,1) + Coorneu(:,2).*Coorneu(:,2));
LR = log(RR);
Theta = atan2(Coorneu(:,2),Coorneu(:,1)); %work in polar coordinates
Phi2theta  = zeros(Nbpt,1);

 J0 = find(Refneu==1);
 J1 = find(Refneu==3);
 J2 = find(Refneu==4);%metal
 J3 = find(Refneu==5);
 
 J4 = find(Refneu==6);
 J5 = find(Refneu==7);%vaccum
 J6 = find(Refneu==8);
 J7 = find(Refneu==2);

 JM = unique([J1; J2; J3; J0]); %dof metal
 JV = unique([J4; J5; J6; J7]); %dof vacuum

 figure;
 subplot(1,2,1)
 scatter(RR(JV).*cos(Theta(JV)),RR(JV).*sin(Theta(JV)),'b');
 hold on;
 scatter(RR(JM).*cos(Theta(JM)),RR(JM).*sin(Theta(JM)),'r');
 hold on;
 axis('equal')
 subplot(1,2,2)
 scatter(LR(JV),Theta(JV),'b');
 hold on;
 scatter(LR(JM),Theta(JM),'r');
 hold on;

switch contrast
    case 'ke < -1'
        
        for k = 1: Nbpt
            
            if (ismember(k,JV) == 1)
                
                if (Theta(k) > 0)
                    Phi2theta(k) = sinh(RacineR*(pi-Theta(k)))/ sinh(RacineR*(pi-phi/2));
                    
                else
                    Phi2theta(k) = -sinh(RacineR*(pi+Theta(k))) / sinh(RacineR*(pi-phi/2));
                end
            end
            if (ismember(k,JM) == 1)
                Phi2theta(k) = sinh(RacineR*(Theta(k)) ) / sinh(RacineR*phi/2.);
            end
            
        end
        
        
    case 'ke > -1'
        
        for k = 1: Nbpt
            
            if (ismember(k,JV) == 1)
                
                if (Theta(k) > 0)
                    Phi2theta(k) = cosh(RacineR*(pi-Theta(k)))/ cosh(RacineR*(pi-phi/2));
                    
                else
                    Phi2theta(k) = cosh(RacineR*(pi+Theta(k))) / cosh(RacineR*(pi-phi/2));
                end
            end
            if (ismember(k,JM) == 1)
                Phi2theta(k) = cosh(RacineR*(Theta(k)) ) / cosh(RacineR*phi/2.);
            end
            
        end        
end

Phi2theta =  RR.^(1i*RacineR).*Phi2theta;
figure;
subplot(1,2,1)
trisurf(Numtri2,LR,Theta,real(Phi2theta));
shading interp;
view(2);
grid off;
colorbar
title(['Singularity in the strip'],'Interpreter','Latex')
subplot(1,2,2)
trisurf(Numtri2,Coorneu(:,1),Coorneu(:,2),real(Phi2theta));
shading interp;
view(2);
grid off;
colorbar
axis('equal')
title(['Singularity at corner'],'Interpreter','Latex')

%%---------------Calcul du coefficient de singularité ---------------------
%Val = 2 * RacineR .* dot(Phi2theta,SS*Phi2theta);
% int 1/e phi^2 Valeur à comparer avec celle obtenue avec Maple
%b = - k^2 * dot(ones(Nbpt,1), SSm * ZZ)/val;
end
%%
% %         Ctrans = (cosh(RacineR*(pi-phi/2))/cosh(-RacineR*phi/2.));
% %         for k = 1: Nbpt
% %            
% %            if (ismember(k,JV) == 1) && (Coorneu(k,1) >=0) 
% %                Phi2theta(k) = cosh(RacineR*(Theta(k)-pi/2));
% %            end
% %            if (ismember(k,JV) == 1) && (Coorneu(k,1) <0) && (Coorneu(k,2) > 0)
% %                Phi2theta(k) = cosh(RacineR*(Theta(k)-pi/2));
% %            end
% %            if (ismember(k,JV) == 1) && (Coorneu(k,1) <0) && (Coorneu(k,2) <= 0)
% %                Phi2theta(k) = cosh(RacineR*(Theta(k)+3*pi/2));
% %            end
% % 
% %         end       
% %         Phi2theta(JM) = Ctrans*cosh(RacineR*(Theta(JM) +(pi/2)));


