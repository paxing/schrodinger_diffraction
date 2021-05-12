% *****************************************************************
% auteurs: Paul Xing, Félix Desrochers, Raphael Thibault
% numéro d'équipe: 10
% date: hiver 2019
% Cours: PHS3903
% *****************************************************************



% **************************Description****************************
% Ce programme permet de réaliser une étude de la conservation de la norme
% du paquet d'onde pour le schéma de différences finies de Crank-Nicolson
% *****************************************************************



%% initialisation des paramètres et du domaine de simulation

clear all, clc


% paramètres ajustables
dt=0.01;
dx=0.05;
dy=dx;
L=15;
H=10;
x_c=L/3;
y_c=H/2;
w=1;
kx=2*sqrt(11);
ky=0;


%initialisation
Nx=round(L/dx)+1;
Ny=round(H/dy)+1;
x=linspace(0,L,Nx);
y=linspace(0,H,Ny);
N=Nx*Ny;
psi1=wave_packet2d(x,y,x_c, y_c,w,w,kx,ky);
psi1vec=reshape(psi1,[N 1]);
psi0vec=psi1vec;

%potentiel nul pour le premier cas
potentiel=sparse(Ny,Nx);

%positionnement de la barrière de potentiel pour le deuxième cas
xp1=ceil(2*Nx/3)-4;
xp2=ceil(2*Nx/3)+4; 
potentiel2=sparse(Ny,Nx);
%barrière infinie modélisée comme étant une fonctionne échelon de très
%haute amplitude 10^10
potentiel2(:,xp1:xp2)=10^10*ones(Ny,xp2-xp1+1);

% calcul matrice pour propagation libre et barrière de potentiel
[M1,M2]=crank_nicolson_R_2d(dt,dx,Nx,Ny,potentiel);
[M1P,M2P]=crank_nicolson_R_2d(dt,dx,Nx,Ny,potentiel2);

%% Enregistrement des résultats en fichier vidéo


psi1vec=psi0vec;
psi1pvec=psi0vec;
vec_t=[];
norme_psi=[];
norme_psi_p=[];


%initialisation du fichier video
x0=10;
y0=10;
width=1350;
height=700;
h = figure('position',[x0,y0,width,height],'visible', 'off');
axis tight manual 
filename = 'analyse_norme_k=6,63_w=1_dx=0,05_dt=0,01.mp4';
v = VideoWriter(filename,'MPEG-4');
v.FrameRate = 8;
v.Quality = 100;
open(v);


for t=1:110
    vec_t=[vec_t t];

    
    subplot(6,2,[1 3 5])
    psi2vec=M2\(M1*psi1vec);
    psi2=reshape(psi2vec, [Ny Nx]);
    prob_num=psi2.*conj(psi2);
    ff=pcolor(x,y,prob_num);
    set(ff, 'EdgeColor', 'none');
    colormap jet;
    c1=colorbar;
    c1.TickLabelInterpreter = 'latex';
    c1.LineWidth = 1;
    title('$$ |\psi_{num}|^2 $$ avec propagation libre', 'Interpreter','latex')
    xlabel('x $$  (a_0) $$','Interpreter','latex')
    ylabel('y $$  (a_0) $$','Interpreter','latex')
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    hold off
    
    
    subplot(6,2,[2 4 6]) 
    psi2pvec=M2P\(M1P*psi1pvec);
    psi2p=reshape(psi2pvec, [Ny Nx]);
    prob_num_p=psi2p.*conj(psi2p);
    ff=pcolor(x,y,prob_num_p);
    set(ff, 'EdgeColor', 'none');
    colormap jet;
    c1=colorbar;
    c1.TickLabelInterpreter = 'latex';
    c1.LineWidth = 1;
    title('$$ |\psi_{num}|^2 $$ avec barri{\`e}re de potentiel', 'Interpreter','latex')
    xlabel('x $$  (a_0) $$','Interpreter','latex')
    ylabel('y $$  (a_0) $$','Interpreter','latex')
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    line([xp2*dx xp2*dx], [10 0],'Color',[1 1 1],'LineWidth',1.5);
    line([xp1*dx xp1*dx], [10 0],'Color',[1 1 1],'LineWidth',1.5);
    hold off
    

    subplot(6,2,[9 11])   
    norme_psi=[norme_psi abs(1-trapz(x,trapz(y,prob_num)))];
    semilogy(vec_t,norme_psi,'Color',[0 0.25 0.75],'Linewidth',1.5)
    xlabel('Nombre de pas de temps','Interpreter','latex');
    ylabel('Erreur','Interpreter','latex')
    ylim([1e-17 1e-14]);
    title('Erreur sur la norme= $$\left|1-\int_{}\int |\psi_{num}|^2dxdy\right|$$','Interpreter','latex')
    legend({strcat('Temps =  ', num2str(dt*t,'%0.2f'),' u.a.')},'Interpreter','latex')
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
   
    
    
    subplot(6,2,[10 12])
    norme_psi_p=[norme_psi_p abs(1-trapz(x,trapz(y,prob_num_p)))];
    semilogy(vec_t,norme_psi_p,'Color',[0 0.25 0.75],'Linewidth',1.5)
    xlabel('Nombre de pas de temps','Interpreter','latex');
    ylabel('Erreur','Interpreter','latex')
    ylim([1e-17 1e-14]);
    title('Erreur sur la norme= $$\left|1-\int_{}\int |\psi_{num}|^2dxdy\right|$$','Interpreter','latex')
    legend({strcat('Temps =  ', num2str(dt*t,'%0.2f'),' u.a.')},'Interpreter','latex')
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    

      
      
    frame = getframe(h);
    writeVideo(v,frame);
    
    
    psi1vec=psi2vec;
    psi1pvec=psi2pvec;
    t=t
end
close(gcf)
close(v);

%% Enregistrement des résultats finaux sous forme de figure

h=figure();
plot(vec_t(1:100),norme_psi(1:100),'r','Linewidth',4);
hold on
plot(vec_t(1:100),norme_psi_p(1:100),'Color',[0 0.25 0.75],'Linewidth',4)
xlabel('Nombre de pas de temps','Interpreter','latex');
ylabel('Erreur sur la norme','Interpreter','latex')
line([50 50], [1e-17 1e-15],'Color','k','LineStyle','--','LineWidth',4);
ylim([1e-17 1e-15]);
%title('Erreur sur la norme= $$\left|1-\int_{}\int |\psi_{num}|^2dxdy\right|$$','Interpreter','latex')
legend({'Propagation libre', 'Avec barri\`{e}re de potentiel','Contact avec la barri\`{e}re'},'Interpreter','latex')
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',32)
set(gca,'FontName','cmr12')
x0=10;
y0=10;
width=1200;
height=600;
set(gcf,'position',[x0,y0,width,height])

print(h,'erreur_norme_v2','-dpng','-r300')
