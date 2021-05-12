% *****************************************************************
% auteurs: Paul Xing, Félix Desrochers, Raphael Thibault
% numéro d'équipe: 10
% date: hiver 2019
% Cours: PHS3903
% *****************************************************************

% **************************Description****************************
% ce programme contient le code principal pour la génération d'une
% animation de la diffraction d'un paquet d'onde gaussien à travers deux
% fentes
% *****************************************************************


%% Initialisation

clear all, clc

% paramètres ajustables du paquet d'onde
dt=0.01;
dx=0.05;
dy=dx;
L=25;
H=25;
x_c=5;
y_c=0;
wx=2.5;
wy=4;
kx=20;
ky=0;


%initialisation du domaine de simulation
Nx=round(L/dx)+1;
Ny=round(H/dy)+1;
x=linspace(0,L,Nx);
y=linspace(-H/2,H/2,Ny);
N=Nx*Ny;

%initialisation du paquet d'onde
psi1=wave_packet2d(x,y,x_c, y_c,wx,wy,kx,ky);
psi1vec=reshape(psi1,[N 1]);
psi0vec=psi1vec;


% définition du potentiel des fentes
potentiel=sparse(Ny,Nx);
% position de la double fente
PF=10; %distance en u.a
epsilon=2; %épaisseur des fentes
xp1=ceil(PF/dx)-epsilon/2;
xp2=ceil(PF/dx)+epsilon/2; 


% paramètre des fentes (taille en nombre de pas dx)

% NOTE: pour modéliser la diffraction à travers une seule fente, il suffit
% de mettre le paramètre a=0.
a=15; %a=0 pour diffraction à travers une fente

b=25;
yp1=ceil(Ny/2)+a;
yp2=ceil(Ny/2)+b;
yp3=ceil(Ny/2)-b;
yp4=ceil(Ny/2)-a;
% ligne de projection en 1D selon y
line_1=20; %distance en u.a.

% la barrière est modélisée comme un une fonction échelon de grande
% amplitude 10^10 (puits infini)
potentiel(:,xp1:xp2)=10^10*ones(Ny,xp2-xp1+1);
potentiel(yp1:yp2,xp1:xp2)=0;
potentiel(yp3:yp4,xp1:xp2)=0;

%appel de la fonction pour initialiser les matrices des coefficients
[M1,M2]=crank_nicolson_R_2d(dt,dx,Nx,Ny,potentiel);

%% Animation de la diffraction du paquet d'onde

%initialisation pour le premier pas de temps
psi1vec=psi0vec;
psi1pvec=psi0vec;
vec_t=[];

%initialisation du fichier video
x0=10;
y0=10;
width=1600;
height=700;
h = figure('position',[x0,y0,width,height],'visible', 'off');
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'modelisation_2_fentes_a=0,75_b=1,25_e=0,1_25x25_k=20_wx=2,5_wy=2,5_dx=0,05_dt=0,01.mp4';
v = VideoWriter(filename,'MPEG-4');
v.FrameRate = 6;
v.Quality = 100;
open(v);



T=170;%nombre d'itération de la simulation
Mat_projection1=zeros(Ny,T); %initalisation de la matrice pour enregister
% les projections sur y


% NOTE: mettre en commentaires les lignes mentionnées dans la boucle pour
% la diffraction à travers 1 fente.
for t=1:T
    %mise à jour du vecteur temps
    vec_t=[vec_t t];
    
    psi2vec=M2\(M1*psi1vec);
    psi2=reshape(psi2vec, [Ny Nx]);
    prob_num=psi2.*conj(psi2);
    
    
    %paquet d'onde en 2D
    subplot(2,6,[1 2 3 7 8 9])
    f1=pcolor(x,y,prob_num);
    set(f1, 'EdgeColor', 'none');
    colormap jet;
    c1=colorbar;
    c1.TickLabelInterpreter = 'latex';
    c1.LineWidth = 1;
    xlabel('x $$  (a_0) $$','Interpreter','latex')
    ylabel('y $$  (a_0) $$','Interpreter','latex')
    dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color',...
        'none');
    legend(dummyh,{strcat('Temps =  ',' ', num2str(dt*t,'%0.2f'),' u.a.')},...
        'Interpreter','latex','TextColor','w','Color','none', ...
        'EdgeColor','none') 
    set(gca,'FontSize',21)
    set(gca,'TickLabelInterpreter','latex')
    %title('Profil de diffraction', 'Interpreter','latex')
    ax=gca;
    ax.LineWidth=1.5;

   
    %dessin des fentes
    
    line([xp1*dx xp1*dx], [H/2 yp2*dx-H/2],'Color',[1 1 1],'LineWidth',1.5);

    line([xp1*dx xp1*dx], [yp3*dx-H/2 -H/2],'Color',[1 1 1],'LineWidth',1.5);
     
    line([xp2*dx xp2*dx], [H/2 yp2*dx-H/2],'Color',[1 1 1],'LineWidth',1.5);
    
    line([xp2*dx xp2*dx], [yp3*dx-H/2 -H/2],'Color',[1 1 1],'LineWidth',1.5);
 
    line([xp1*dx xp2*dx], [yp2*dx-H/2 yp2*dx-H/2],'Color',[1 1 1],...
        'LineWidth',1.5);
    
    line([xp1*dx xp2*dx], [yp3*dx-H/2 yp3*dx-H/2],'Color',[1 1 1],...
        'LineWidth',1.5);
   
     
    
    
    % NOTE: pour modéliser la diffraction à travers une seule fente, mettre
    % ces lignes en commentaires
    % ******************************************************
    line([xp1*dx xp1*dx], [yp1*dx-H/2 yp4*dx-H/2],'Color',[1 1 1],...
       'LineWidth',1.5);
    line([xp2*dx xp2*dx], [yp1*dx-H/2 yp4*dx-H/2],'Color',[1 1 1],...
       'LineWidth',1.5);
    line([xp1*dx xp2*dx], [yp1*dx-H/2 yp1*dx-H/2],'Color',[1 1 1],...
       'LineWidth',1.5);
    line([xp1*dx xp2*dx], [yp4*dx-H/2 yp4*dx-H/2],'Color',[1 1 1],...
       'LineWidth',1.5);
    % ******************************************************

    
    
    % dessin de la ligne de projection
    line([line_1 line_1],[-H/2 H/2],'Color','red','LineStyle','--',...
        'LineWidth',1.5)
    projection_1=prob_num(:,ceil(line_1/dx));
    Mat_projection1(:,t)=projection_1;

    
    %projection sur y
    subplot(2,6,[ 11 12])
    plot(y,projection_1,'Color',[0 0.25 0.75],'Linewidth',2)
    xlabel('y $$  (a_0) $$','Interpreter','latex')
    ylabel(' $$  |\psi_{num}|^2 $$','Interpreter','latex')
    title(' Projection 1D apr{\`e}s les 2 fentes', 'Interpreter','latex')
    set(gca,'FontSize',21)
    set(gca,'TickLabelInterpreter','latex')
    ylim([0 0.25e-1])
    xlim([-H/2 H/2])
    ax=gca;
    ax.LineWidth=1.5;
    
    % densité de probabilité en 3d
    subplot(2,6,[5 6])
    s=surf(x,y,prob_num);
    s.EdgeColor = 'none';
    colormap jet;
    c2=colorbar;
    c2.TickLabelInterpreter = 'latex';
    c2.LineWidth = 1;
    zlim([0 0.25e-1])
    xlim([PF L]);
    ylim([-H/2 H/2])
    zlabel(' $$  |\psi_{num}|^2 $$','Interpreter','latex')
    caxis([0 0.25e-1])
    title(' $$ |\psi_{num}|^2 $$ apr{\`e}s les 2 fentes', 'Interpreter','latex')
    ax=gca;
    ax.LineWidth=1.5;
    set(gca,'FontSize',21)
    set(gca,'TickLabelInterpreter','latex')

    
    set(gcf,'position',[x0,y0,width,height])
    frame = getframe(h);
    writeVideo(v,frame);
    psi1vec=psi2vec;
    
    t=t
end
close(gcf)
close(v);


%% 
% ************************************************************
% NOTE: LES SECTIONS SUIVANTES NE S'APPLIQUE QUE POUR LA DIFFRACTION À 
% TRAVERS DEUX FENTES
% ************************************************************

%% régression avec le modèle théorique de diffraction à travers 2 fentes

% patron final de diffraction selon y de la simulation
Diffraction_sim=Mat_projection1(:,end);
Diffraction_sim=Diffraction_sim-min(Diffraction_sim);

%régression avec paramettre a,b,c
%note: il peut être nécessaire d'ajuster les «start point» au cas par cas
fitresult = fit(y',Diffraction_sim,'c*(sinc((b-a)*x))^2*(cos((a+b)*x))^2',...
    'StartPoint', [a*dx b*dx max(Diffraction_sim)]);
coef= coeffvalues(fitresult);
coef_int=confint(cfit(fitresult));
coef_err=(coef_int(2,:)-coef_int(1,:))/2;

%% Figure de la comparaison du patron de diffraction analytique et numérique

h=figure(2)
plot(y,Diffraction_sim,'Color',[0 0.25 0.75],'Linewidth',2)
hold on

fit1=coef(3)*(sinc((coef(2)-coef(1))*y)).^2.*(cos((coef(1)+coef(2))*y)).^2;
plot(y,fit1,'r','Linewidth',2)
legend({'Simulation','Th\''{e}orique'},'Interpreter','latex')
xlabel('y $$  (a_0) $$','Interpreter','latex')
ylabel(' $$  |\psi_{num}|^2 $$','Interpreter','latex')
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',28)
set(gca,'FontName','cmr12')
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
print(h,'fit_double_fente','-dpng','-r0')


%% Figure avec représentation 2d de la simulation et de la comparaison


h=figure()


%représentation 2d
subplot(2,7,[1 2 3 8 9 10])
f1=pcolor(x,y,prob_num);
set(f1, 'EdgeColor', 'none');
colormap jet;
c1=colorbar;
c1.TickLabelInterpreter = 'latex';
c1.LineWidth = 1;
caxis([0 1.1*max(Diffraction_sim)])
zlim([0 1.1*max(Diffraction_sim)])
xlim([5 25])
xlabel('x $$  (a_0) $$','Interpreter','latex')
ylabel('y $$  (a_0) $$','Interpreter','latex')
set(gca,'FontSize',28)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontName','cmr12')
ax=gca;
ax.LineWidth=2;

%dessin des fentes

line([xp1*dx xp1*dx], [H/2 yp2*dx-H/2],'Color',[1 1 1],'LineWidth',1.5);

line([xp1*dx xp1*dx], [yp3*dx-H/2 -H/2],'Color',[1 1 1],'LineWidth',1.5);

line([xp2*dx xp2*dx], [H/2 yp2*dx-H/2],'Color',[1 1 1],'LineWidth',1.5);

line([xp2*dx xp2*dx], [yp3*dx-H/2 -H/2],'Color',[1 1 1],'LineWidth',1.5);

line([xp1*dx xp2*dx], [yp2*dx-H/2 yp2*dx-H/2],'Color',[1 1 1],...
    'LineWidth',1.5);

line([xp1*dx xp2*dx], [yp3*dx-H/2 yp3*dx-H/2],'Color',[1 1 1],...
    'LineWidth',1.5);





% ******************************************************
line([xp1*dx xp1*dx], [yp1*dx-H/2 yp4*dx-H/2],'Color',[1 1 1],...
    'LineWidth',1.5);
line([xp2*dx xp2*dx], [yp1*dx-H/2 yp4*dx-H/2],'Color',[1 1 1],...
    'LineWidth',1.5);
line([xp1*dx xp2*dx], [yp1*dx-H/2 yp1*dx-H/2],'Color',[1 1 1],...
    'LineWidth',1.5);
line([xp1*dx xp2*dx], [yp4*dx-H/2 yp4*dx-H/2],'Color',[1 1 1],...
    'LineWidth',1.5);
% ******************************************************


%dessin de la ligne de projection
line([line_1 line_1],[-H/2 H/2],'Color','red','LineStyle','--','LineWidth',1.5)

%comparaison de la projection sur y avec le modèle analytique
subplot(2,7,[5 6 7 12 13 14])
plot(y,Diffraction_sim,'Color',[0 0.25 0.75],'Linewidth',2)
hold on
fit=coef(3)*(sinc((coef(2)-coef(1))*y)).^2.*(cos((coef(1)+coef(2))*y)).^2;
plot(y,fit,'r','Linewidth',2)
legend({'Simulation','Mod\`{e}le th\''{e}orique'},'Interpreter','latex')
xlabel('y $$  (a_0) $$','Interpreter','latex')
ylabel(' $$  |\psi_{num}|^2 $$','Interpreter','latex')
ylim([0 1.1*max(Diffraction_sim)])
ax=gca;
xlim([-12.5 12.5])
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',28)
set(gca,'FontName','cmr12')
x0=10;
y0=10;
width=2100;
height=700;
set(gcf,'position',[x0,y0,width,height])

print(h,'fit_double_fente_test2','-dpng','-r0')

%% Analyse de l'effet de la variation des paramètres a et b


clear all, clc

% variables fixés
dt=0.01; dx=0.05; dy=dx;
L=25; H=25;
x_c=5; y_c=0;
wx=2.5; wy=6;
kx=20; ky=0;
Nx=round(L/dx)+1; Ny=round(H/dy)+1;
x=linspace(0,L,Nx); y=linspace(-H/2,H/2,Ny);
N=Nx*Ny;
psi1=wave_packet2d(x,y,x_c, y_c,wx,wy,kx,ky);
psi1vec=reshape(psi1,[N 1]);
psi0vec=psi1vec;
potentiel=sparse(Ny,Nx);
PF=10; epsilon=2;
xp1=ceil(PF/dx)-epsilon/2; xp2=ceil(PF/dx)+epsilon/2; 
line_1=20;

% paramètre des fentes 
a_vec=[10 15 15 15 15 17 17];
b_vec=[25 20 25 30 35 32 35];

% fonction de régression: 'c*(sinc((b-a)*x))^2*(cos((a+b)*x))^2'
fit_param_coef=zeros(3,length(a_vec));
% avec comme paramètre respectif:
%fit_param_coef(1,:)=a
%fit_param_coef(1,:)=b
%fit_param_coef(1,:)=c

fit_param_error=zeros(3,length(a_vec));
Diffraction_mat=zeros(length(y),length(a_vec));
% seule la projection sur y à la fin de la simulation est enregistrée
for i=1:length(a_vec)
    a=a_vec(i);
    b=b_vec(i);
    
    yp1=ceil(Ny/2)+a;
    yp2=ceil(Ny/2)+b;
    yp3=ceil(Ny/2)-b;
    yp4=ceil(Ny/2)-a;
    potentiel(:,xp1:xp2)=10^10*ones(Ny,xp2-xp1+1);
    potentiel(yp1:yp2,xp1:xp2)=0;
    potentiel(yp3:yp4,xp1:xp2)=0;
    [M1,M2]=crank_nicolson_R_2d(dt,dx,Nx,Ny,potentiel);
    
    psi1vec=psi0vec;
    psi1pvec=psi0vec;
    vec_t=[];
    T=170;%nombre d'itération de la simulation
 
    
    for t=1:T
        psi2vec=M2\(M1*psi1vec);
        psi1vec=psi2vec;
    end
    psi2=reshape(psi2vec, [Ny Nx]);
    prob_num=psi2.*conj(psi2);
    Diffraction_sim=prob_num(:,ceil(line_1/dx))
    Diffraction_sim=Diffraction_sim-min(Diffraction_sim);
    Diffraction_mat(:,i)=Diffraction_sim;
    
    
    %régression des paramètres
    
    fitresult = fit(y',Diffraction_sim,'c*(sinc((b-a)*x))^2*(cos((a+b)*x))^2',...
    'StartPoint', [a*dx/2 b*dx/2 max(Diffraction_sim)]);
    coef= coeffvalues(fitresult);
    coef_int=confint(cfit(fitresult));
    coef_err=(coef_int(2,:)-coef_int(1,:))/2;
    fit_param_coef(:,i)=coef';
    fit_param_error(:,i)=coef_err';
end

%note: il peut être nécessaire d'ajuster les «start point» au cas par cas
% auquel cas, la matrice  Diffraction_mat permet de récuper le profil de
% diffraction pour chaque cas (b,a). Il est alors possible d'ajuster
% manuellement les régressions si nécessaire.
%% Figure de la relation b'/a' vs b/a



a_fit=fit_param_coef(1);
b_fit=fit_param_coef(2);
a_err=fit_param_error(1,:);
b_err=fit_param_error(2,:);



ratio1=b_vec./a_vec;
ratio2=b_fit./a_fit;
err_ratio=ratio2.*(a_err./a_fit+b_err./b_fit);



h=figure()

%régression linéaire y=mx+b de la relation entre les paramètres
coef_pente=polyfit(ratio1,ratio2,1);
fit1=coef_pente(2)+coef_pente(1).*linspace(1.2,2.8,31);
plot(linspace(1.2,2.8,31),fit1,'r','Linewidth',2)
hold on
scatter(ratio1,ratio2,100,'b','filled')
hold on
errorbar(ratio1,ratio2,err_ratio, 'LineStyle','none','LineWidth',2, 'Color','b')
ax=gca;
ax.LineWidth=2;

xlabel('$$b/a$$','Interpreter','latex')
ylabel('$$b''/a''$$','Interpreter','latex')
xlim([1 3])
legend1=strcat('$$y =$$', num2str(coef_pente(1),'%0.4f'),'$$x + $$',num2str(coef_pente(2),'%0.3f'));
legend({legend1,'Donn\''{e}es'},'Interpreter','latex','Location','NW')
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',28)
set(gca,'FontName','cmr12')
x0=10;
y0=10;
width=900;
height=600;
set(gcf,'position',[x0,y0,width,height])
print(h,'fit_double_fente_parameter_fit_v2','-dpng','-r0')


%%





