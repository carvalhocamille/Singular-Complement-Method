function [M1,M2] =NCQ(Coorneu,Numtri,Reftri,Gd,G,Gds,Gs,mur_metal,mur_dielec,kw)
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

for t=1:nt %loop over triangles
    
    %Elementary matrices
    S=[Coorneu(Numtri(t,1),:);Coorneu(Numtri(t,2),:);Coorneu(Numtri(t,3),:)];
    S21=S(2,:)-S(1,:);S31=S(3,:)-S(1,:);
    delta=S21(1)*S31(2)-S21(2)*S31(1);
    
    %Elementary matrices by quadratures
    Mt1=zeros(q);Mt2=zeros(q);
    
    for k=1:nbq %loop over the triangle DOFs 
        x=pts_quadT(k,1);y=pts_quadT(k,2);
        
        if(q==3)  %P1
            w=[1-x-y x y];
        else      % P2
            w=[(1-x-y)*(1-2*x-2*y) x*(2*x-1) y*(2*y-1) 4*x*(1-x-y) 4*x*y 4*y*(1-x-y)];
        end
        
        pk=weights_closed(k)*abs(delta); % weight by delta in physical domain
        
        if (Reftri(t)==3 || Reftri(t)==4  || Reftri(t)==5)
            Mt1= Mt1 +  pk * w' * ( Gd(Numtri(t,k)) + kw^2 * mur_metal * G(Numtri(t,k)) );  
            Mt2= Mt2 +  pk * w' * ( Gds(Numtri(t,k)) + kw^2 * mur_metal * Gs(Numtri(t,k)) );
        else
            Mt1= Mt1 + pk * w' * ( Gd(Numtri(t,k)) + kw^2 * mur_dielec * G(Numtri(t,k)) );
            Mt2= Mt2 + pk * w' * ( Gds(Numtri(t,k))+ kw^2 * mur_dielec * Gs(Numtri(t,k)) );
            
        end
        
    end
    
    for i=1:q
        I=Numtri(t,i);
        M1(I)=M1(I)+Mt1(i);
        M2(I)=M2(I)+Mt2(i);
    end
    
end  
end

