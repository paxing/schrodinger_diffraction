% *****************************************************************
% auteurs: Paul Xing, F�lix Desrochers, Raphael Thibault
% num�ro d'�quipe: 10
% date: hiver 2019
% Cours: PHS3903
% *****************************************************************



% **************************Description****************************
% Cette fonction permet de calculer la propagation th�orique de la densit�
% de probabilit� d'un paquet d'onde gaussien sym�trique en x et y.
% *****************************************************************


% **************************Param�tres******************************
% Entr�e:
%   x = vecteur d�finissant le maillage en x
%   y = vecteur d�finissant le maillage en y
%   x_c = position initiale du centre du paquet d'onde en x
%   y_x = position initiale du centre du paquet d'onde en y
%   w = largeur du paquet d'onde (w=2 delta x)
%   kx = nombre d'onde en x 
%   ky = nombre d'onde en y
%   t = temps de propagation 

% Sortie:
%   prob: matrixe y par x de la densit� de probabilit� du paquet d'onde
% *****************************************************************


function prob= wave_packet_prob_th(x,y,x_c, y_c,w,kx,ky,t)

%initialisation de la matrice d�finissant la densit� de probabilit�
M=length(x);
N=length(y);
prob=zeros(N,M);

% �cart-type du paquet d'onde d�pendant du temps
sigma_t2=w^2/4+t^2/w^2;


for i=1:M
    for j=1:N
        prob(j,i)=1/(2*pi*sigma_t2)*exp(-(x(i)-x_c-kx*t)^2/(2*sigma_t2))...
            *exp(-(y(j)-y_c-ky*t)^2/(2*sigma_t2));
    end
end

