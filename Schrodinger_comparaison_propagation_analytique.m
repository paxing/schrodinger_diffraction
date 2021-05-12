% *****************************************************************
% auteurs: Paul Xing, Félix Desrochers, Raphael Thibault
% numéro d'équipe: 10
% date: hiver 2019
% Cours: PHS3903
% *****************************************************************


% **************************Description****************************
% Ce programme permet de générer une animation comparant la propagation
% d'un paquet d'onde à l'aide des différences finies (Crank-Nicolson) avec
% la propagation analytique du paquet d'onde
% *****************************************************************



%% initialisation

clear all, clc

dt=0.01;
dx=0.05;
dy=dx;
L=15;
H=7.5;
x_c=L/3;
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

%% enregistrement video de la comparaison

psi1vec=psi0vec;
psi1pvec=psi0vec;
vec_t=[];
norme_psi=[];
norme_psi_p=[];


%initialisation du fichier vidéo
x0=10;
y0=10;
width=900;
height=1400;
h = figure('position',[x0,y0,width,height],'visible', 'off');
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'analyse_comparaison_th_k=6,63_w=1_dx=0,05_dt=0,01.mp4';
v = VideoWriter(filename,'MPEG-4');
v.FrameRate = 6;
v.Quality = 100;
open(v);


for t=1:100
    vec_t=[vec_t t];

    
    subplot(3,1,[1])
    
    %modèle mumérique
    psi2vec=M2\(M1*psi1vec);
    psi2=reshape(psi2vec, [Ny Nx]);
    prob_num=psi2.*conj(psi2);
    f1=pcolor(x,y,prob_num);
    set(f1, 'EdgeColor', 'none');
    colormap jet;
    c1=colorbar;
    c1.TickLabelInterpreter = 'latex';
    c1.LineWidth = 1;
    title('Mod{\`e}le num{\''e}rique $$ |\psi_{num}|^2 $$', 'Interpreter','latex')
    xlabel('x $$  (a_0) $$','Interpreter','latex')
    ylabel('y $$  (a_0) $$','Interpreter','latex')
    set(gca,'FontSize',16)
    set(gca,'TickLabelInterpreter','latex')
    hold off
    
    
    
    subplot(3,1,[2])
    
    % modèle théorique
    prob_th= wave_packet_prob_th(x,y,x_c, y_c,w,kx,ky,t*dt);
    f2=pcolor(x,y,abs(prob_th));
    set(f2, 'EdgeColor', 'none');
    colormap jet;
    c2=colorbar;
    c2.TickLabelInterpreter = 'latex';
    c2.LineWidth = 1;
    title('Mod{\`e}le th{\''e}orique $$ |\psi_{th}|^2 $$', 'Interpreter','latex')
    xlabel('x $$(a_0)$$', 'Interpreter','latex')
    ylabel('y $$(a_0)$$',  'Interpreter','latex')
    set(gca,'FontSize',16)
    hold off

 

    subplot(3,1,[3])
    
    %erreur
    f3=pcolor(x,y,abs(prob_th-prob_num));
    set(f3, 'EdgeColor', 'none');
    colormap jet;
    c2=colorbar;
    c2.TickLabelInterpreter = 'latex';
    c2.LineWidth = 1;
    title('$$ \left||\psi_{th}|^2-|\psi_{num}|^2\right| $$', 'Interpreter','latex')
    xlabel('x $$(a_0)$$', 'Interpreter','latex')
    ylabel('y $$(a_0)$$',  'Interpreter','latex')
    set(gca,'FontSize',16)
    hold off

    frame = getframe(h);
    writeVideo(v,frame);
    psi1vec=psi2vec;
    t=t
end
close(gcf)
close(v);
%% Enregistrement d'une video version 2 pour le ppt

psi1vec=psi0vec;
psi1pvec=psi0vec;
vec_t=[];
norme_psi=[];
norme_psi_p=[];

%initialisation du fichier
x0=10;
y0=10;
width=2400;
height=800;
h = figure('position',[x0,y0,width,height],'visible', 'off');
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'analyse_comparaison_v2_th_k=6,63_w=1_dx=0,05_dt=0,01.mp4';
v = VideoWriter(filename,'MPEG-4');
v.FrameRate = 6;
v.Quality = 100;
open(v);


for t=1:100
    vec_t=[vec_t t];

    subplot(2,3,[1])
    
    %modèle mumérique
    psi2vec=M2\(M1*psi1vec);
    psi2=reshape(psi2vec, [Ny Nx]);
    prob_num=psi2.*conj(psi2);
    f1=pcolor(x,y,prob_num);
    set(f1, 'EdgeColor', 'none');
    colormap jet;
    c1=colorbar;
    c1.TickLabelInterpreter = 'latex';
    c1.LineWidth = 1;
    title('Mod{\`e}le num{\''e}rique $$ |\psi_{num}|^2 $$', 'Interpreter','latex')
    xlabel('x $$  (a_0) $$','Interpreter','latex')
    ylabel('y $$  (a_0) $$','Interpreter','latex')
    set(gca,'FontSize',32)
    set(gca,'TickLabelInterpreter','latex')
    hold off
    
    subplot(2,3,[4])
    
    % modèle théorique
    prob_th= wave_packet_prob_th(x,y,x_c, y_c,w,kx,ky,t*dt);
    f2=pcolor(x,y,abs(prob_th));
    set(f2, 'EdgeColor', 'none');
    colormap jet;
    c2=colorbar;
    c2.TickLabelInterpreter = 'latex';
    c2.LineWidth = 1;
    title('Mod{\`e}le th{\''e}orique $$ |\psi_{th}|^2 $$', 'Interpreter','latex')
    xlabel('x $$(a_0)$$', 'Interpreter','latex')
    ylabel('y $$(a_0)$$',  'Interpreter','latex')
    set(gca,'FontSize',32)
    hold off

    %erreur
    pos1 = [0.45 0.25 0.4 0.6];
    subplot(2,3,[2 3 5 6],'Position',pos1)
    f3=pcolor(x,y,abs(prob_th-prob_num));
    set(f3, 'EdgeColor', 'none');
    colormap jet;
    c2=colorbar;
    c2.TickLabelInterpreter = 'latex';
    c2.LineWidth = 1;
    title('$$ \left||\psi_{th}|^2-|\psi_{num}|^2\right| $$', 'Interpreter','latex')
    xlabel('x $$(a_0)$$', 'Interpreter','latex')
    ylabel('y $$(a_0)$$',  'Interpreter','latex')
    set(gca,'FontSize',32)
    hold off
    
    set(gcf,'position',[x0,y0,width,height])

    
    frame = getframe(h);
    writeVideo(v,frame);

    psi1vec=psi2vec;
   t=t
end
close(gcf)
close(v);


%% Sauvagarde sous forme de figure de la propagation finale

h=figure();

%modèle mumérique
subplot(3,1,[1])
psi2vec=M2\(M1*psi1vec);
psi2=reshape(psi2vec, [Ny Nx]);
prob_num=psi2.*conj(psi2);
f1=pcolor(x,y,prob_num);
set(f1, 'EdgeColor', 'none');
colormap jet;
c1=colorbar;
c1.TickLabelInterpreter = 'latex';
c1.LineWidth = 2;
title('Mod{\`e}le num{\''e}rique $$ |\psi_{num}|^2 $$', 'Interpreter','latex')
xlabel('x $$  (a_0) $$','Interpreter','latex')
ylabel('y $$  (a_0) $$','Interpreter','latex')
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',18)
set(gca,'FontName','cmr12')


% modèle théorique
subplot(3,1,[2])
prob_th= wave_packet_prob_th(x,y,x_c, y_c,w,kx,ky,t*dt);
f2=pcolor(x,y,abs(prob_th));
set(f2, 'EdgeColor', 'none');
colormap jet;
c2=colorbar;
c2.TickLabelInterpreter = 'latex';
c2.LineWidth = 2;
title('Mod{\`e}le th{\''e}orique $$ |\psi_{th}|^2 $$', 'Interpreter','latex')
xlabel('x $$(a_0)$$', 'Interpreter','latex')
ylabel('y $$(a_0)$$',  'Interpreter','latex')
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',18)
set(gca,'FontName','cmr12')

%erreur
subplot(3,1,[3])
f3=pcolor(x,y,abs(prob_th-prob_num));
set(f3, 'EdgeColor', 'none');
colormap jet;
c2=colorbar;
c2.TickLabelInterpreter = 'latex';
c2.LineWidth = 2;
title('$$ \left||\psi_{th}|^2-|\psi_{num}|^2\right| $$', 'Interpreter','latex')
xlabel('x $$(a_0)$$', 'Interpreter','latex')
ylabel('y $$(a_0)$$',  'Interpreter','latex')
ax=gca;
ax.LineWidth=2;
set(gca,'TickLabelInterpreter','latex','FontName','cmr12')
set(gca,'FontSize',18)
set(gca,'FontName','cmr12')

x0=10;
y0=10;
width=600;
height=900;
set(gcf,'position',[x0,y0,width,height])
print(h,'camparaison_modele_analytique_v2','-dpng','-r300')

