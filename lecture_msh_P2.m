%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lecmail.m:
% routine de lecture de fichiers de maillages triangulaires 2D au format .msh 
% 
% SYNOPSIS [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecmail(nomfile)
%          
% INPUT * nom_maillage : le nom d'un fichier de maillage au format msh
%                   AVEC SON SUFFIXE .msh
%
% OUTPUT - Nbpt : nbre de sommets (entier) 
%        - Nbtri : nbre de triangles (entier)
%        - Coorneu : coordonnees (x, y) des sommets (matrice reelle Nbpt x 2)
%        - Refneu : reference des sommets (vecteur entier Nbpt x 1)
%        - Numtri : liste de triangles 
%                   (3 numeros de sommets - matrice entiere Nbtri x 6)
%        - Reftri : reference des triangles (matrice entiere Nbtri x 1)
%        - Nbaretes : Nombre d'aretes sur la frontiere.
%        - Numaretes(Nbaretes,3) = numero des 2 Noeuds de chaque Arete de la frontiere 
%		 - Refaretes(Nbaretes,1) = Reference de chaque ar?te de la
%		 fronti?re
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh_P2(nomfile)

fid=fopen(nomfile,'r');
if fid <=0,
   msg=['Le fichier de maillage : ' nomfile ' n''a pas ?t? trouv?'];
   error(msg);
end
while ~strcmp(fgetl(fid),'$Nodes'), end
Nbpt = str2num(fgetl(fid));
Coorneu = zeros(Nbpt,2);
Refneu = zeros(Nbpt,1);
Nbaretes = 0;
Numaretes = [];
Refaretes = [];
RefneuBis = zeros(Nbpt,1);
for i=1:Nbpt
tmp= str2num(fgetl(fid));
Coorneu(i,:) = tmp(2:3);
end
while ~strcmp(fgetl(fid),'$Elements'), end
Nbtri = str2num(fgetl(fid));
tmp= str2num(fgetl(fid)); test = tmp(2);
%while test~=15, tmp= str2num(fgetl(fid)); test = tmp(2); Nbtri = Nbtri-1;end
% on traite les noeuds des coins
k=1;
Nbelem=Nbtri;
while test~=9
    RefneuBis(tmp(end-2:end)')=tmp(4);
    Numaretes= [Numaretes;tmp(end-2:end)];
    Refaretes= [Refaretes;tmp(4)];
    tmp= str2num(fgetl(fid));
    test = tmp(2); Nbaretes=Nbaretes+1; Nbtri=Nbtri-1; k=k+1;
end
Numtri = zeros(Nbtri,6);
Reftri = zeros(Nbtri,1);
while test~=8
    Refneu(tmp(end-5:end)')=tmp(4);    %attention refneu c'est pas reftri
    Numtri(k-Nbaretes,:)=tmp(end-5:end)';
    Reftri(k-Nbaretes)=tmp(4);
    k=k+1;
    if (k<=Nbelem)
        tmp= str2num(fgetl(fid)); test = tmp(2);
    else
        break;
    end
    %Nbtri = Nbtri-1;
end
% on traite les noeuds du bord
Refneu(find(RefneuBis==1))=RefneuBis(find(RefneuBis==1));
Refneu(find(RefneuBis==2))=RefneuBis(find(RefneuBis==2));
fclose(fid);