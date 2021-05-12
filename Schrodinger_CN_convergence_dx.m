% *****************************************************************
% auteurs: Paul Xing, Félix Desrochers, Raphael Thibault
% numéro d'équipe: 10
% date: hiver 2019
% Cours: PHS3903
% *****************************************************************

% **************************Description****************************
% Ce programme permet de réaliser une étude de la convergence en dx pour la
% propagation d'un paquet d'onde utilisant le schéma de différences finies
% de Crank-Nicolson
% *****************************************************************


%% Calcul de l'erreur pour différents pas dx

clear all, clc

% pas de temps fixe
dt=0.005;

%choix des pas en dx
vec_dx=[0.025 0.05 0.1 0.2];


% enregistrement de l'erreur dans des cellules
cell_erreur=cell(1,length(vec_dx));
cell_norme=cell(1,length(vec_dx));

%parcours les différents pas dx 
for i=1:length(vec_dx)
    
    %initialisation du maillage pour chaque dx
    dx=vec_dx(i);
    dy=dx;
    L=20;
    H=20;
    x_c=L/4;
    y_c=H/2;
    w=1;
    kx=2*sqrt(11);
    ky=0;
    
    Nx=round(L/dx)+1;
    Ny=round(H/dy)+1;
    
    x=linspace(0,L,Nx);
    y=linspace(0,H,Ny);
    N=Nx*Ny;
    
    psi1=wave_packet2d(x,y,x_c, y_c,w,w,kx,ky);
    psi1vec=reshape(psi1,[N 1]);
    psi0vec=psi1vec;
    potentiel=sparse(Ny,Nx);
    [M1,M2]=crank_nicolson_R_2d(dt,dx,Nx,Ny,potentiel);
    psi1vec= psi0vec;
    vec_t=[];
    norme_psi=[];
    erreur_psi=[];
    
    % nombre d'itération temporelle fixe pour chaque pas dx
    tn=100;
        for t=1:tn
           
            psi2vec=M2\(M1*psi1vec);
            psi2=reshape(psi2vec, [Ny Nx]);
            prob_num=psi2.*conj(psi2);
            prob_th= wave_packet_prob_th(x,y,x_c, y_c,w,kx,ky,t*dt);
            erreur_psi=[erreur_psi (trapz(x,trapz(y,abs(prob_th-prob_num))))];
            norme_psi=[norme_psi abs(1-trapz(x,trapz(y,prob_num)))];
            psi1vec=psi2vec;
        end
        
        cell_erreur{i}= erreur_psi;
        cell_norme{i}= norme_psi;
        
        disp(i)
end

%% figure convergence dx

 % cette section permet de tracer l'évolution de l'erreur en fonction du
 % nombre de pas de temps pour les différents pas d'espace choisis

h=figure();
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
 for i=1:length(vec_dx)
    v_t=dt:dt:dt*tn; 
    semilogy(v_t, cell_erreur{i},'-o','Linewidth',2.5,'MarkerSize',8,...
        'MarkerFaceColor','auto');
    hold on
 end
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',32)
set(gca,'FontName','cmr12')
xlabel('Temps (u.a.)','Interpreter','latex','FontName','cmr12')
ylabel('Erreur absolue','Interpreter','latex','FontName','cmr12')
title('Erreur =$$ \int_{}\int\left||\psi_{th}|^2-|\psi_{num}|^2\right|dxdy $$',...
    'Interpreter','latex')
legend({'$$dx=0.025$$','$$dx=0.05$$', '$$dx=0.1$$', '$$dx=0.2$$'},...
    'Location','SE','Interpreter','latex','FontName','cmr12')



print(h,'CN_convergence_dx_dt=0,005_v3','-dpng','-r300')

%% ordre de convergence dx

% cette section permet de calculer l'ordre de convergence et de tracer la
% figure de la relation log(erreur) vs log(dx)

% permet de prendre l'erreur finale pour chaque pas d'espace dx
err_f=[];
for i=1:length(vec_dx)
    err_f(i)=cell_erreur{i}(end);
end

%permet d'obtenir les coefficients de regression dont la pente détermine
%l'ordre de convergence
coef=coeffvalues(fit(vec_dx',err_f','power1'));

h=figure();
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
loglog(vec_dx(1:end),err_f,'b-o','Linewidth',4,'MarkerSize',15,...
    'MarkerFaceColor','auto')
hold on
xlabel('$$dx (a_0)$$','Interpreter','latex','FontName','cmr12')
ylabel('Erreur absolue','Interpreter','latex','FontName','cmr12')
xlim([vec_dx(1) vec_dx(end)])


f1 = coef(1)*vec_dx.^(coef(2));
plot(vec_dx,f1,'r-','Linewidth',4)
hold on
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',32)
set(gca,'FontName','cmr12')
legend({ strcat('m =  ',' ', num2str(coef(2),'%0.3f')),'R{\''e}gression'},...
    'Interpreter','latex')
print(h,'Ordre_CN_convergence_dx_dt=0,005_v2','-dpng','-r0')