% *****************************************************************
% auteurs: Paul Xing, Félix Desrochers, Raphael Thibault
% numéro d'équipe: 10
% date: hiver 2019
% Cours: PHS3903
% *****************************************************************

% **************************Description****************************
% Ce programme permet de réaliser une étude de la convergence en dt pour la
% propagation d'un paquet d'onde utilisant le schéma de différences finies
% de Crank-Nicolson
% *****************************************************************


%% initialisation du maillage et des paramètres du paquet d'onde


clear all, clc

% paramètres ajustables
dx=0.025;
dy=dx;
L=20;
H=20;
x_c=L/4;
y_c=H/2;
w=1;
kx=2*sqrt(11);
ky=0;


% initialisation du domaine de simulation
Nx=round(L/dx)+1;
Ny=round(H/dy)+1;
x=linspace(0,L,Nx);
y=linspace(0,H,Ny);
N=Nx*Ny;
psi1=wave_packet2d(x,y,x_c, y_c,w,w,kx,ky);
psi1vec=reshape(psi1,[N 1]);
psi0vec=psi1vec;
potentiel=sparse(Ny,Nx);


%% Calcul de l'erreur

% variation du pas de temps pour le calcul de convergence
 dt= [0.005 0.01 0.02 0.04];
 
 % le temps total de la simulation a été fixé à 1 u.a.
 
 % calcul du nombre de pas de temps T (attention aux choix de pas de temps
 % car T doit être entier)
 T=1*ones(1,length(dt))./dt;
 
 %initialisation des cellules pour enregistrer l'erreur des différents pas
 %de temps
 cell_erreur=cell(1,length(dt));
 cell_norme=cell(1,length(dt));
 
 %parcours les différents pas de temps choisi
 for i=1:length(dt)
     
     % initialisation des matrices pour chaque pas de temps
     [M1,M2]=crank_nicolson_R_2d(dt(i),dx,Nx,Ny,potentiel);
     psi1vec=psi0vec;
     vec_t=[];
     norme_psi=[];
     erreur_psi=[];
     
     %calcul pour T itérations
     for t=1:T(i)
         psi2vec=M2\(M1*psi1vec);
         psi2=reshape(psi2vec, [Ny Nx]);
         prob_num=psi2.*conj(psi2);
         prob_th= wave_packet_prob_th(x,y,x_c, y_c,w,kx,ky,t*dt(i));
         erreur_psi=[erreur_psi (trapz(x,trapz(y,abs(prob_th-prob_num))))];
         norme_psi=[norme_psi abs(1-trapz(x,trapz(y,prob_num)))];
         psi1vec=psi2vec;
     end
      cell_erreur{i}= erreur_psi;
      cell_norme{i}= norme_psi;
      disp(i)
 end
 

 %% figure convergence dt

 % cette section permet de tracer l'évolution de l'erreur en fonction du
 % nombre de pas de temps pour les différents pas de temps choisi
 
 
 h=figure(1)
 
 x0=10;
 y0=10;
 width=900;
 height=600;
 set(gcf,'position',[x0,y0,width,height])
 
 for i=1:length(dt)
    v_t=dt(i):dt(i):1; 
    semilogy(v_t, cell_erreur{i},'-o','Linewidth',2.5,'MarkerSize',8,...
        'MarkerFaceColor','auto');
    hold on
 end
 
 ax=gca;
 ax.LineWidth=2;
 set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
 set(gca,'FontSize',32)
 set(gca,'FontName','cmr12')
 xlabel('Temps (u.a.)','Interpreter','latex','FontName','cmr12')
 ylabel('Erreur absolue','Interpreter','latex','FontName','cmr12')
 title('Erreur =$$ \int_{}\int\left||\psi_{th}|^2-|\psi_{num}|^2\right|dxdy $$',...
     'Interpreter','latex')
 legend({'$$dt = 0.005 $$', '$$dt = 0.01$$', '$$dt = 0.02$$','$$dt = 0.04 $$'},...
     'Location','SE','Interpreter','latex')
 print(h,'CN_convergence_dt_dx=0,025_v2','-dpng','-r300')


%% convergence dt


% cette section permet de calculer l'ordre de convergence et de tracer la
% figure de la relation log(erreur) vs log(dt)

err_f=[];
norm_f=[];

% permet de prendre l'erreur finale pour chaque pas de temps dt
for i=1:length(dt)
    err_f(i)=cell_erreur{i}(end);
    norm_f(i)=cell_norme{i}(end);
end

%permet d'obtenir les coefficients de regression dont la pente détermine
%l'ordre de convergence

coef=coeffvalues(fit(dt',err_f','power1'));

h=figure();
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
loglog(dt,err_f,'b-o','Linewidth',4,'MarkerSize',15,'MarkerFaceColor','auto')
hold on
f1 = coef(1)*dt.^(coef(2));
plot(dt,f1,'r-','Linewidth',4)
hold on
xlabel('$$dt$$ (u.a.)','Interpreter','latex')
ylabel('Erreur absolue','Interpreter','latex','FontName','cmr12')
xlim([dt(1) dt(end)])
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',32)
set(gca,'FontName','cmr12')
legend({ strcat('m =  ',' ', num2str(coef(2),'%0.3f')),'R{\''e}gression'},'Interpreter','latex')
print(h,'Ordre_CN_convergence_dt_dx=0,025_v2','-dpng','-r0')