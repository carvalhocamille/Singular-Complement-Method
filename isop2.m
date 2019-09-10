%
%           [Num1]=isop2(Num2)
% 
% conversion d'un maillage P2 en maillage P1 (iso)
% d�coupage d'un triangle en 4 triangles suivant les points milieux des ar�tes
% 
%  Num2 : tableau nt x 6 de la num�rotation P2
%
%  Num1 : tableau 4nt x 3 de la num�rotation P1

function [Num1]=isop2(Num2)
nt=size(Num2,1);
Num1=zeros(4*nt,3);
k=0;
for t=1:nt,
    k=k+1; Num1(k,1)=Num2(t,1);Num1(k,2)=Num2(t,4);Num1(k,3)=Num2(t,6);
    k=k+1; Num1(k,1)=Num2(t,4);Num1(k,2)=Num2(t,2);Num1(k,3)=Num2(t,5);
    k=k+1; Num1(k,1)=Num2(t,4);Num1(k,2)=Num2(t,5);Num1(k,3)=Num2(t,6);
    k=k+1; Num1(k,1)=Num2(t,3);Num1(k,2)=Num2(t,6);Num1(k,3)=Num2(t,5);
end;
    