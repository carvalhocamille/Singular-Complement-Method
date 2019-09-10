function [M1,M2] =NCQ(Coorneu,Numtri,Reftri,Gd,G,Gds,Gs,mur_metal,mur_dielec,epsilonr_metal, epsilonr_dielec)
% Compute the Newton-Cotes Quadrature (NCQ) over triangles.
% The weights and nodes for second order NCQ can be found in "Symmetric
% Quadrature Formulae for Simplexes", P. Silvester


%Formula with 6 points, second order, the points being the DOF for
%Lagrange P2

%Nodes and weight on the reference triangle
pts_quadT=[0 0; 1 0; 1 0; 0.5 0; 0.5 0.5; 0 0.5];
weights_closed = [0; 0; 0; 1/3 ; 1/3; 1/3];
weights_open = [7/12; 7/12; 7/12; -1/4 ; -1/4; -1/4];
  
%Number of triangles and DOFs
nt=size(Numtri,1);
ns=size(Coorneu,1);
q=size(Numtri,2);
nbq=length(weights_open);

M1=zeros(ns,1);
M2=zeros(ns,1);

for t=1:nt
    
    %Elementary matrices
    S=[Coorneu(Numtri(t,1),:);Coorneu(Numtri(t,2),:);Coorneu(Numtri(t,3),:)];
    S21=S(2,:)-S(1,:);S31=S(3,:)-S(1,:);
    delta=S21(1)*S31(2)-S21(2)*S31(1);
    Jflmt=[S31(2) -S21(2);-S31(1) S21(1)]/delta; %transfo affine!
    
    %Elementary matrices by quadratures
    Mt1=zeros(q);Mt2=zeros(q); Mt11=zeros(q);Mt22=zeros(q);
    for k=1:nbq
        x=pts_quadT(k,1);y=pts_quadT(k,2);
        
        if(q==3)  %P1
            w=[1-x-y x y];
            gw=[-1 1 0;-1 0 1];
        else % P2
            w=[(1-x-y)*(1-2*x-2*y) x*(2*x-1) y*(2*y-1) 4*x*(1-x-y) 4*x*y 4*y*(1-x-y)]; 
            gw=[4*(x+y)-3 4*x-1 0 4*(1-2*x-y) 4*y -4*y;
                4*(x+y)-3 0 4*y-1 -4*x 4*x 4*(1-x-2*y)];
        end
 
        %coordonn�es dans le plan physique
        P=S'*[1-x-y; x; y];

        pk=weights_closed(k)*abs(delta);
        
        if Reftri(t)==3 || Reftri(t)==4  || Reftri(t)==5 
            Mt1=Mt1 +  pk * w' * Gd(k);
            Mt11=Mt11 +  mur_metal * pk * w' * G(k);
            Mt2=Mt2 +  pk * w' * Gds(k);
            Mt22=Mt22 +  mur_metal * pk * w' * Gs(k);
           % jg=Jflmt*gw;
           % Kt=Kt+(1./epsilonr_metal) * pk * jg' * Gd(k); %gives a vector of size (6,1)
        else
            Mt1=Mt1 + pk * w' * Gd(k);
            Mt11=Mt11 + mur_dielec * pk * w' * G(k);
            Mt2=Mt2 + pk * w' * Gds(k);
            Mt22=Mt22 + mur_dielec * pk * w' * Gs(k);
           % jg=Jflmt*gw;
           % Kt=Kt+(1./epsilonr_dielec) * pk * jg' * Gd(k);
        end

    end
    
    %assemblage de M et K
   % I = Numtri(t,:)';
   % M1(I) = M1(I)+ Mt1;
   % M2(I) = M2(I) + Mt2;
    
    for i=1:q
        I=Numtri(t,i);
        M1(I)=M1(I)+Mt1(i) + Mt11(i);
        M2(I)=M2(I)+Mt2(i) + Mt22(i);
    end
    
end  %fin boucle t
end

