% *****************************************************************
% auteurs: Paul Xing, Félix Desrochers, Raphael Thibault
% numéro d'équipe: 10
% date: hiver 2019
% Cours: PHS3903
% *****************************************************************

% **************************Description****************************
% Cette fonction permet d'initialiser un paquet d'onde gaussien
% *****************************************************************


% **************************Paramètres******************************
% Entrée:
%   x = vecteur définissant le maillage en x
%   y = vecteur définissant le maillage en y
%   x_c = position initiale du centre du paquet d'onde en x
%   y_x = position initiale du centre du paquet d'onde en y
%   wx = largeur du paquet d'onde en x(w=2 delta x)
%   wy = largeur du paquet d'onde en y(w=2 delta y)
%   kx = nombre d'onde en x 
%   ky = nombre d'onde en y


% Sortie:
%   prob: matrixe y par x définissant la fonction d'onde
% *****************************************************************


function psi= wave_packet2d(x,y,x_c, y_c,wx,wy,kx,ky)


% initialisation de la matrice
M=length(x);
N=length(y);
psi=zeros(N,M);

for i=1:M
    for j=1:N
        % calculer du paquet d'onde comme étant factorisable en x et y
        psi(j,i)=exp(1.0i*kx*x(i))*exp(-(x(i)-x_c)^2/wx^2)*exp(1.0i*ky*y(j))*exp(-(y(j)-y_c)^2/wy^2);
    end
end

% calcul de la constante de normalisation (permet d'initialiser le paquet 
% d'onde avec une fonction d'onde normalisée)
f=conj(psi).*psi;
Int=trapz(x,trapz(y,f));
c_norm=1/sqrt(Int);

psi=c_norm*psi;