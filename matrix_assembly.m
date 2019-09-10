%
%         [K,M]=calcul_EF_2D(Corneu,Numtri,Reftri,fK,fM)
%
%calcul des matrices de rigidit� et de masse EF P1 ou P2 (�l�ment droit)
%
%  Corneu : matrice ns x 2 des coordonn�es des noeuds du maillage 
%  Numtri : matrice nt x nd de la num�rotation globale des triangles
%          (ns=nb de noeuds, nt=nb de trianngles, nd=nb noeuds par triangle)
%  Reftri : vecteur nt x 1 des r�f�rences par triangle
%  fK,fM  : fonctions scalaires de (x,y,ref) associ�es aux termes 
%           de rigidit� et de masse (prise en compte de mat�riaux diff�rents)
%
%  K  : matrice sparse nsxns de rigidit� volumique
%  M  : matrice sparse nsxns de masse volumique

function [KK,MM,SS,SSm]=matrix_assembly(Coorneu,Numtri,Reftri,Numaretes,Refarete,Nbaretes,epsilonr_dielec,epsilonr_metal,mur_dielec,mur_metal,kw)

%formule de quadrature 2D sur le triangle de r�f�rence
%-----------------------------------------------------
%Formule de Hammer-Stroud a 3 points,exacte pour P2 sur le triangle (STROUD p. 307)
  pts_quadT=[0 0.5;0.5 0;0.5 .5];
  os=1/6;pds_quadT=os*[1 1 1];
%Formule de Radon, Hammer, Marlowe & Stroud a 7 points, exacte pour P5 sur le triangle (STROUD p.314)
% os=sqrt(15);s3=1./3.;
% pp1=(6.-os)/21.;pp2=(6.+os)/21.;pp3=(9.+2.*os)/21.;pp4=(9.-2.*os)/21.;
% pts_quadT=[s3 s3;pp1 pp1;pp1 pp3;pp3 pp1;pp2 pp2;pp2 pp4;pp4 pp2];
% pp1=(155.-os)/2400.;pp2=(155.+os)/2400.;
% pds_quadT=[9./80.;pp1;pp1;pp1;pp2;pp2;pp2];

%tailles diverses
nt=size(Numtri,1);
ns=size(Coorneu,1);
q=size(Numtri,2);  % 3 en P1 ou 6 en P2
nbq=length(pds_quadT);
%formule de quadrature 1D, ordre 3
s35=sqrt(3./5);
pts_quadS=0.5*[1-s35 1 1+s35];
osS=1/18;pds_quadS=osS*[5 8 5];
nbqS=length(pds_quadS);

%d�claration des matrices
KK=sparse(ns,ns);
MM=sparse(ns,ns);
SS=sparse(ns,ns);
SSm=sparse(ns,ns);

%BOUCLE PRINCIPALE SUR LES TRIANGLES
%-----------------------------------
for t=1:nt
    
    %calcul des matrices �l�mentaires
    S=[Coorneu(Numtri(t,1),:);Coorneu(Numtri(t,2),:);Coorneu(Numtri(t,3),:)];
    S21=S(2,:)-S(1,:);S31=S(3,:)-S(1,:);
    delta=S21(1)*S31(2)-S21(2)*S31(1);
    Jflmt=[S31(2) -S21(2);-S31(1) S21(1)]/delta; %transfo affine!
    
    %calcul des matrices �l�mentaires Mt et Kt par quadrature
    Mt=zeros(q,q);Kt=zeros(q,q);
    for k=1:nbq
        x=pts_quadT(k,1);y=pts_quadT(k,2);
        %calcul des fonctions de base aux points de quadrature
        if(q==3)  %P1
            w=[1-x-y x y];
        %fprintf('taille de w\n');
        %size(w)            
            gw=[-1 1 0;-1 0 1];
        else      %P2
            w=[(1-x-y)*(1-2*x-2*y) x*(2*x-1) y*(2*y-1) 4*x*(1-x-y) 4*x*y 4*y*(1-x-y)];
        % fprintf('taille de w\n');
        %size(w)           
            gw=[4*(x+y)-3 4*x-1 0 4*(1-2*x-y) 4*y -4*y;
                4*(x+y)-3 0 4*y-1 -4*x 4*x 4*(1-x-2*y)];
        end
        %coordonn�es dans le plan physique
        %P=S'*[1-x-y; x; y];
        pk=pds_quadT(k)*abs(delta);
        
        if Reftri(t)==3 || Reftri(t)==4  || Reftri(t)==5 
            Mt=Mt- kw^2 * mur_metal *pk*(w'*w);
            jg=Jflmt*gw;
            Kt=Kt+(1./epsilonr_metal)*pk*(jg'*jg);
        else
            Mt=Mt- kw^2 * mur_dielec *pk*(w'*w);
            jg=Jflmt*gw;
            Kt=Kt+(1./epsilonr_dielec)*pk*(jg'*jg);
        end

    end
    
    %assemblage de M et K
    for i=1:q
        I=Numtri(t,i);
        for j=1:q
            J=Numtri(t,j);
            KK(I,J)=KK(I,J)+Kt(i,j);
            MM(I,J)=MM(I,J)+Mt(i,j);
        end
    end
    
end  %fin boucle t

%BOUCLE PRINCIPALE SUR LES TRIANGLES
for t=1:Nbaretes
  
    %boucle sur les ar�tes
  % for a=1:3,
 
        %as=mod(a,3)+1;
        %bs=mod(a+1,3)+1;
        I=Numaretes(t,1);J=Numaretes(t,2);K=Numaretes(t,3);
        %ar�te sur la fronti�re ref_bord
        %if(ismember(Bord(t,a),ref_bord))       
           % L1=norm(Coorneu(I,:)-Coorneu(J,:)); %longueur ar�te1
            L=norm(Coorneu(I,:)-Coorneu(J,:)); %longueur ar�te2
            %L=L1+L2;
            for k=1:nbq                   % boucle quadrature 1D
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
                end %fin test degr� EF
                
                if Refarete(t)==1
                    St = St *(1./epsilonr_metal);
                    Stm = St * mur_metal;
                else
                    St = St *(1./epsilonr_dielec);
                    Stm = St * mur_dielec;
                end
                
                %assemblage de M
                p=size(in,2);
                for i=1:p
                    for j=1:p
                        SS(in(i),in(j))=SS(in(i),in(j))+St(i,j);
                        SSm(in(i),in(j))=SSm(in(i),in(j))+Stm(i,j);
                    end
                end
            end       
        %end, %fin test frontiere 
    %end,  %fin boucle a
end  %fin boucle t









    