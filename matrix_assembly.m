function [KK,MM,SS,SSm]=matrix_assembly(Coorneu,Numtri,Reftri,Numaretes,Refarete,Nbaretes,epsilonr_dielec,epsilonr_metal,mur_dielec,mur_metal)
% Compute P2 FE matrices using 3points Hammer-Stroud quadrature (STROUD p. 307)

%Formula with 6 points, second order, the points being the DOF for
%Lagrange P2

%Nodes and weight on the reference triangle
pts_quadT=[0 0.5;0.5 0;0.5 .5];
os=1/6;pds_quadT=os*[1 1 1];

nt=size(Numtri,1);
ns=size(Coorneu,1);
q=size(Numtri,2);  
nbq=length(pds_quadT);

%1D quadrature rule, order 3
s35=sqrt(3./5);
pts_quadS=0.5*[1-s35 1 1+s35];
osS=1/18;pds_quadS=osS*[5 8 5];
nbqS=length(pds_quadS);

KK=sparse(ns,ns);
MM=sparse(ns,ns);
SS=sparse(ns,ns);
SSm=sparse(ns,ns);


for t=1:nt %loop over triangles
    
    %Elementary matrices
    S=[Coorneu(Numtri(t,1),:);Coorneu(Numtri(t,2),:);Coorneu(Numtri(t,3),:)];
    S21=S(2,:)-S(1,:);S31=S(3,:)-S(1,:);
    delta=S21(1)*S31(2)-S21(2)*S31(1);
    Jflmt=[S31(2) -S21(2);-S31(1) S21(1)]/delta; % affine transformation 
    
    %Elementary matrices by quadratures
    Mt=zeros(q,q);Kt=zeros(q,q);
    
    for k=1:nbq %loop over the triangle DOFs 
        x=pts_quadT(k,1);y=pts_quadT(k,2);
        
        %basis functions
        if(q==3)          %P1
            w=[1-x-y x y];         
            gw=[-1 1 0;-1 0 1];
        else              %P2
            w=[(1-x-y)*(1-2*x-2*y) x*(2*x-1) y*(2*y-1) 4*x*(1-x-y) 4*x*y 4*y*(1-x-y)];       
            gw=[4*(x+y)-3 4*x-1 0 4*(1-2*x-y) 4*y -4*y;
                4*(x+y)-3 0 4*y-1 -4*x 4*x 4*(1-x-2*y)];
        end

        pk=pds_quadT(k)*abs(delta);
        
        if Reftri(t)==3 || Reftri(t)==4  || Reftri(t)==5 
            Mt=Mt + mur_metal *pk*(w'*w);
            jg=Jflmt*gw;
            Kt=Kt+(1./epsilonr_metal)*pk*(jg'*jg);
        else
            Mt=Mt + mur_dielec *pk*(w'*w);
            jg=Jflmt*gw;
            Kt=Kt+(1./epsilonr_dielec)*pk*(jg'*jg);
        end

    end
    
    %assembly
    for i=1:q
        I=Numtri(t,i);
        for j=1:q
            J=Numtri(t,j);
            KK(I,J)=KK(I,J)+Kt(i,j);
            MM(I,J)=MM(I,J)+Mt(i,j);
        end
    end
    
end  

%For the surface matrices
for t=1:Nbaretes
    
    I=Numaretes(t,1);J=Numaretes(t,2);K=Numaretes(t,3);
    L=norm(Coorneu(I,:)-Coorneu(J,:)); %length
    
    for k=1:nbq                   %  quadrature 1D
        x=pts_quadS(k);
        c=L*pds_quadS(k);
        if(q==3)                     %P1
            w=[1-x x];
            St=c*w'*w;
            in=[I J];
        else                         %P2
            w=[((1-2*x)*(1-x)) ((2*x-1)*x) (4*x*(1-x))];
            St=c*w'*w;
            in=[I J K];
        end
        
        if Refarete(t)==1
            St = St *(1./epsilonr_metal);
            Stm = St * mur_metal;
        else
            St = St *(1./epsilonr_dielec);
            Stm = St * mur_dielec;
        end
        
        %assembly
        p=size(in,2);
        for i=1:p
            for j=1:p
                SS(in(i),in(j))=SS(in(i),in(j))+St(i,j);
                SSm(in(i),in(j))=SSm(in(i),in(j))+Stm(i,j);
            end
        end
    end
end









    