% AE_RES Test Case Ic
% lets give the whole thing a better physical connection by introducing the
% proper order of wind speed.
% T = 1
% dt = 0.001

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%   INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length of domain
L_dom = 100;% [m]
% relativel length of supply zone
alpha = 0.2;

% Threshold velocity
U_th = 4; %[m/s]
% Simulation time
total_time = 1000; %[s]
% adaptation timescale
T=0.01; %[s]
% constant wind
u_const = 7.5;% [m/s]
% u_w(1:total_time,1:101) = 7.5; %[m/s]
VS=1;
z=[];

% source magnitude(s)

src = [0.5 0.7  2 4]*1.5e-4*(u_const-U_th)^3/(L_dom*alpha); %[Kg/m2s]


% numerical stuff
dx=1;dt=T;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_line = 0:dx:L_dom;

source=zeros(total_time/dt,L_dom/dx+1);
u_w(1:total_time/dt,1:L_dom/dx+1) = u_const; 

% Ct (=concentration in transport)
Ct = zeros(total_time/dt,L_dom/dx+1);
% Ca (=concentration available on bed)
Ca = zeros(total_time/dt,L_dom/dx+1);

% Cu =(concentration based on wind speed)
Cu = zeros(total_time/dt,L_dom/dx+1);

figure(23);close;figure(23)
% figure
hold all

plot_style= (['k-.';'k: ';'k--';'k- '])
for i=1:length(src)
    
source(:,2:L_dom*alpha/dx)= src(i)*dt*dx; % wind is in equilibrium with supply if is around 0.01-0.02

    
% [Ct,u_w,Ca,Cu,Ccap_index] = core_upwind_matrix_T(u_w,U_th,source,dx,dt,total_time,T);
[Ct,Ca,Cu,Ccap_index] = model_core(u_w,U_th,source,dx,dt,total_time,T,VS,z);


plot(x_line,Ct(end,:)./Cu(end,:),plot_style(i,:),'linewidth',2)
% plot(x_line,Ct(end,:)./Cu(end,:),'linewidth',2)

%  plot_results(Ct,Ca,u_w,source)

end
legend(sprintf('\\beta = %3.2f',src(1)*L_dom*alpha/(1.5e-4*(u_const-U_th)^3)),...
    sprintf('\\beta = %3.2f',src(2)*L_dom*alpha/(1.5e-4*(u_const-U_th)^3)),...
    sprintf('\\beta = %3.2f',src(3)*L_dom*alpha/(1.5e-4*(u_const-U_th)^3)),...
     sprintf('\\beta = %3.2f',src(4)*L_dom*alpha/(1.5e-4*(u_const-U_th)^3)))

% legend('supply << demand','supply < demand','supply > demand')

% legend('S_s = 0.01 kg/m^{2}s','S_s = 0.02 kg/m^{2}s','S_s = 0.05 kg/m^{2}s')
legend('Location','NEO')
plot([20 20],[0 1.1],'color',[0.8 0.8 0.8])
xlabel('Distance - x [m]')
ylabel('Q/Q_{u} [-]')

xlim([0 110])
ylim([0 1.1])

set(gcf,'PaperUnits','centimeters','PaperPosition',[0.634517 6.34517 12 5]) 
 