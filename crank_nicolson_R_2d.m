% *****************************************************************
% auteurs: Paul Xing, Félix Desrochers, Raphael Thibault
% numéro d'équipe: 10
% date: hiver 2019
% Cours: PHS3903
% *****************************************************************

% **************************Description****************************
% Cette fonction permet de calculer les matrices M1 et M2 des différences
% finies pour l'équation de Schrodinger en 2d dépendante du temps. Elle
% utile le schéma implicite de Crank-Nicolson et la méthode de la matrice
% 2D.
%
% Les conditions frontières utilisées sont de types réféchisantes.
% *****************************************************************


% **************************Paramètres******************************
% Entrée:
%   dt= pas de temps
%   dx = pas d'espace
%   Nx = taille du maillage en x
%   Ny = taille du maillage en y
%   potentiel = matrixe Ny par Nx décrivant le potentiel présent

% Sortie:
%   M1 = matrice des coefficients multipliant psi_n
%   M2 = matrice des coefficients multipliant ps_(n+1)
% *****************************************************************




function [M1,M2]=crank_nicolson_R_2d(dt,dx,Nx,Ny,potentiel)



alpha=1.0i*dt/(4*dx^2);
N=Nx*Ny;

%initialisation des matrices sparses
M1=sparse(N,N);
M2=sparse(N,N);

%initialisation du vecteur v commun à placer sur la diagonale +1 et -1 
%(pour plus de détail, voir annexe du rapport final)

v=[];
%calcul des coefficients du vecteur v
for i=1:Nx-2
    vec=[alpha*ones(1,Ny-2) 0 0];
v =[v vec];
end

%vecteur à placer sur la diagonale -1
v1= [zeros(1, Ny) v zeros(1,Ny)];

%vecteur à placer sur la diagonale +1
v2= [zeros(1, Ny+2) v zeros(1,Ny)];

%ajout des vecteurs dans les matrices des coefficients
M1=spdiags(v1(:),-1,N,N)+spdiags(v2(:),1,N,N);
M2=spdiags(-v1(:),-1,N,N)+spdiags(-v2(:),1,N,N);




% vecteur à placer sur la diagonale -Ny
w1=[0 v zeros(1,Ny)];
% vecteur à placer sur la diagonale +Ny
w2=[zeros(1,2*Ny+1) v];

%ajout des vecteurs dans les matrices des coefficients
M1=M1+spdiags(w1(:),-Ny,N,N)+spdiags(w2(:),Ny,N,N);
M2=M2+spdiags(-w1(:),-Ny,N,N)+spdiags(-w2(:),Ny,N,N);


%transformation du potentiel est des coeff de la diagonale principale en
%vecteur
diag_principale_1=1-1.0i*dt.*reshape(potentiel,[N,1])/2-4*alpha;
diag_principale_2=1+1.0i*dt.*reshape(potentiel,[N,1])/2+4*alpha;
diag_principale_1(1:Ny)=0;
diag_principale_2(1:Ny)=-1;
diag_principale_1(N-Ny:N)=0;
diag_principale_2(N-Ny:N)=-1;

% Ajustement pour tenir compte des conditions frontières réféchisantes
for i=1:Nx-1
    diag_principale_1(Ny+(i-1)*Ny+1)=0;
    diag_principale_1(Ny+i*Ny)=0;
    
    diag_principale_2(Ny+(i-1)*Ny+1)=-1;
    diag_principale_2(Ny+i*Ny)=-1;
end


%ajout des vecteurs dans les matrices des coefficients
M2=M2+spdiags(diag_principale_2(:),0,N,N);
M1=M1+spdiags(diag_principale_1(:),0,N,N);




