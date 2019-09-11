%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SCM.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM to solve scattering problem using P2 Lagrange FE, and with Singular
% Complement Method (SCM)
% The problem solved is the following:
% |  div(1/eps Grad u)+ mu u = 0,        in D,
% |                     u = f,                 on B, 
% where D is a disk composed of 2 materials characterized by
% eps and mu: these parameters can change sign

% This code uses a .msh mesh generated by GMSH that dictates the geoemtry
% The code is organized as follows:
% 1) Define parameters 
% 2) Read mesh. We use lecture_msh_P2.m to extract information from the .msh
% 3) Assemble FE matrices. We use matrix_assembly.m to assemble mass, stiffness, and surface matrices
% 4) Build right hand side and solve with Standard FEM
% 5) Build the SCM: solve the linear system to compute the singular complement 
%       We use cutoff.m to define the cutoff function,
%              singular_expression.m to compute the analytic terms involving the singularity
%              NCQ.m to compute the coupling terms via Newton-Cotes
%              Quadrature rule
% 6) Compute the singularity coefficient (b)
% 7) Solve the problem using SCM
% 8) save and post-process the results

%More information can be found in the paper, 
% "The Singular Complement Method for metamaterial-dielectric transmission
% problems" Carvalho, Ciarlet (2019).
% Camille Carvalho, PhD
% 9/10/2019
%%
 clear all;
 close all;
 
% =====================================================
%% 1) Parameters 
% =====================================================

x0 = 0;                              %center of domain (x0,y0)
y0 = 0;
epsilonP = 1;                        %P for Positive material              
epsilonM = -0.83;                    %M for Metamaterial
muP = 1;                   
muM = 1; 
angle= pi/6;                         %angle of the domain
%We chose for the singularity:
RacineR = 1.322404244;               %Singular exponant value 
kw = 3;                              %wavenumber
kappa_epsilon  = epsilonM/ epsilonP; %ratio of permittivities
kappa_mu  = muM/ muP;                %ration of permeabilities


fprintf('Parameter : kappa_epsilon =  %d \n',kappa_epsilon);

% =====================================================
%% 2) Read mesh
% =====================================================
mesh = 'mesh/disk.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes] = lecture_msh_P2(mesh);
Numtri2 = isop2(Numtri);

% Display mesh and test zones
Index0 = find(Refneu==1);
Index1 = find(Refneu==2);
Index2 = find(Refneu==3);
Index3 = find(Refneu==4);
Index4 = find(Refneu==5);
Index5 = find(Refneu==6);
Index6 = find(Refneu==7);
Index7 = find(Refneu==8);

figure(1);
scatter(Coorneu(Index1,1),Coorneu(Index1,2));
hold on;
scatter(Coorneu(Index2,1),Coorneu(Index2,2));
hold on;
scatter(Coorneu(Index3,1),Coorneu(Index3,2));
hold on;
scatter(Coorneu(Index4,1),Coorneu(Index4,2));
hold on;
scatter(Coorneu(Index5,1),Coorneu(Index5,2));
hold on;
scatter(Coorneu(Index6,1),Coorneu(Index6,2));
hold on;
scatter(Coorneu(Index7,1),Coorneu(Index7,2));
axis 'equal'
title(['Mesh and sub domains'],'Interpreter','Latex')

fprintf('Display mesh and zones ... \n')
% =====================================================
%% 3) Matrices assembly
% =====================================================
MM    = sparse(Nbpt,Nbpt);  % mass matrix
KK    = sparse(Nbpt,Nbpt);  % stiffness matrix
SS    = sparse(Nbpt,Nbpt);  % surface mass matrix
F     = ones(Nbpt,1);       % rhs
FFreg = zeros(Nbpt,1);      % rhs for regular part
UU    = zeros(Nbpt,1);      % unknown using std FEM 
UUreg = zeros(Nbpt,1);      % unknown using std FEM + SCM
ZZ    = zeros(Nbpt,1);      % dual singularity (of the form z_h + c_h s_b + s)
X0 = x0*ones(Nbpt,1);
Y0 = y0*ones(Nbpt,1);

%Compute stiffness, mass, and surface FE matrices
[KK,MM,SS,SSm]=matrix_assembly(Coorneu,Numtri,Reftri,Numaretes,Refaretes,Nbaretes,epsilonP,epsilonM,muP,muM,kw);
% =====================================================
%% ----------------Right-hand side---------------------
% =====================================================
FF = -MM*F;

% =====================================================
%% 4) Solving problem with standard FEM
% =====================================================
%Elimination of the Dirichlet boundary condition
bdy = unique([Numaretes(:,1); Numaretes(:,2); Numaretes(:,3)]);
interior = find(ismember([1:Nbpt],bdy)==0);
NbptD = Nbpt - length(bdy);

UUD = zeros(NbptD,1);
KKD = KK(interior,interior);
MMD = MM(interior,interior);
FFD = FF(interior);
fprintf('Start solving for standard FEM.....\n');
AAD =KKD+MMD; 
UUD = AAD\FFD;
UU(interior) = UUD;
fprintf('Solved\n');

% =====================================================
%% 5) Build the SCM 
% =====================================================
% a)Create cutoff function
delta = 0.1;
l= 0.9;
h = 0.001;
phi = pi/6;
[eta, etap, etapp] = cutoff(delta, l, h);

% b)Create exact expression for terms involving singularities
[Beta_s,div_s_grad_s, s_grad_s,Zsingu,sPhi,Phi] = singular_expression(eta, etap, etapp, RacineR,phi,1./epsilonM,1./epsilonP);
[Beta_s_b,div_s_grad_s_b, s_grad_s_b,Zsingu_b,sPhib,Phib] = singular_expression(eta, etap, etapp, -RacineR,phi,1./epsilonM,1./epsilonP);

% c)Create singularity vectors associated with DOFs
div_s_grad_sh = zeros(Nbpt,1);
Zsinguh = zeros(Nbpt,1);
div_s_grad_sh_b = zeros(Nbpt,1);
Zsinguh_b = zeros(Nbpt,1);
for k = 1: Nbpt
    dist = sqrt(Coorneu(k,1)^2 + Coorneu(k,2)^2);
    angle = atan2(Coorneu(k,2),Coorneu(k,1));
    Zsinguh(k) =  Zsingu(dist,angle);
    Zsinguh_b(k) = Zsingu_b(dist,angle);
     if (dist < delta)
        div_s_grad_sh(k) = 0;                       %Impose 0 near the corner
        div_s_grad_sh_b(k) = 0;
    else
        div_s_grad_sh(k)=  div_s_grad_s(dist,angle);
        div_s_grad_sh_b(k)= div_s_grad_s_b(dist,angle);
    end
end

figure;
 trisurf(Numtri2,Coorneu(:,1),Coorneu(:,2),real(Zsinguh));
 colorbar;
 view(2);shading interp;grid off; axis 'equal';
 title('\eta s_h on the mesh')

figure;
 trisurf(Numtri2,Coorneu(:,1),Coorneu(:,2),real(div_s_grad_sh));
 colorbar;
 view(2);shading interp;grid off; axis 'equal'; title('div(\sigma \nabla s)_h on the mesh')
%% d) First approach to compute the dual singularity
AA1 = sparse(NbptD,NbptD);
FF1 = zeros(NbptD,1);
ZZ1 = zeros(NbptD,1);

%Create regular terms with only FEM, need to find the corner
or = find(Coorneu ==[0 0]);
origin = or(1);
inter  = find(ismember(interior,origin)==0);
interior1 = interior(inter);

%Create coupling terms with NCQ
[BdS,BdSb] = NCQ(Coorneu,Numtri,Reftri,div_s_grad_sh,Zsinguh_b,div_s_grad_sh_b,Zsinguh_b,muM,muP,kw);

%Create system to solve to find the dual singularity
AA1(1:end,1:end-1)= -KK(interior,interior1)+MM(interior,interior1);    %Matrix for regular part of the dual singularity, z_h
AA1(:,end) = BdS(interior);                                            %Part including the dual singularity coefficient c_h
FF1 = -BdSb(interior);                                                 % Right-hand side (see of (17) in CC19)

figure; spy(AA1);
figure;
trisurf(Numtri2,Coorneu(:,1),Coorneu(:,2),real(BdSb));
colorbar;
view(2);shading interp;grid off; axis 'equal'; title('real part of right hand side')

%% e) Compute the singularity coefficient
fprintf('Start solving for dual singularity.....\n'); 
ZZ1 = AA1\FF1;
ZZ = zeros(Nbpt,1);
ZZ(interior1) = ZZ1(1:end-1);
Dual_singularity =  ZZ1(end)* Zsinguh_b +Zsinguh_b + ZZ ; 
fprintf('Solved\n');

figure;
trisurf(Numtri2,Coorneu(:,1),Coorneu(:,2),real(Dual_singularity));
colorbar;
shading interp;grid off; axis 'equal'; colorbar; 
 title('dual singularity')
% =====================================================
%% 6) Compute the singularity coefficient (b)
% =====================================================
 [rhs,~] = NCQ(Coorneu,Numtri,Reftri,Dual_singularity,zeros(Nbpt,1),zeros(Nbpt,1),zeros(Nbpt,1),muM,muP,kw);
 Nsteps =100;
 value = 0;
 for k=1:Nsteps
     i = -pi + 2*pi/(Nsteps-1) *(k-1);
     value = value + 2 * 2*pi/(Nsteps-1) * sPhi(i)*Phi(i);
 end
 value = value * 2 * 1i*RacineR;
 b = -dot(F,rhs) /value; %check indices
% =====================================================
%% 7) Solve the problem using SCM
% =====================================================
 
 FFreg = FF(interior1) - b * BdSb(interior1); %check indices
 %FFDreg = FFreg(interior1);
 AAD =  -KK(interior1,interior1)+MM(interior1,interior1);
 fprintf('Start solving for SCM + FEM.....\n'); 
 UUreg = AAD\FFreg;
 fprintf('Solved\n');

% =====================================================
%% 8) Saving
% =====================================================
%TBD
% =====================================================
%% 9) Display results
% =====================================================%
%a) Results from standard FEM
UU = UU + F;                  %add lift
figure;
trisurf(Numtri2,Coorneu(:,1),Coorneu(:,2),real(UU));
colorbar;
view(2);shading interp;grid off; axis 'equal';
 %UM = max(real(UU));
 %Um = min(real(UU));
 %caxis([Um UM]);
caxis([-7 2]);
title(['re(u) std FEM'],'Interpreter','Latex')

%a) Results from FEM + SCM
UUscm = ones(Nbpt,1); 
UUscm(interior1) =  UUreg;
figure;
trisurf(Numtri2,Coorneu(:,1),Coorneu(:,2),real(UUscm));
colorbar;
view(2);shading interp;grid off; axis 'equal';
 %UM = max(real(UU));
 %Um = min(real(UU));
 %caxis([Um UM]);
caxis([-7 2]);
title(['re(u) SCM + FEM'],'Interpreter','Latex')


