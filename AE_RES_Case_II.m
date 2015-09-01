% AE_RES Test Case II
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%   INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length of domain
L_dom = 100;% [m]
% relativel length of supply zone
alpha = 0.2;
% source magnitude(s)
src = [0.45e-3]; %[Kg/m2s]
% Threshold velocity
U_th = 4; %[m/s]
% Simulation time
total_time = 90; %[s]
% adaptation timescale
T=0.5; %[s]

VS=1;
z=[];


% numerical stuff
dx=1;dt=0.01;


% varying wind
u_w = zeros(total_time/dt,L_dom/dx+1);

u_w(1:40/dt,:)=8;
u_w((40/dt+1):55/dt,:)=6;
u_w((55/dt+1):90/dt,:)=9;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




source=zeros(total_time/dt,L_dom/dx+1);
source(:,2:L_dom*alpha/dx)= src*dt*dx;

% alter wind conditions
if 0; u_w(1:total_time/2,:)=u_w(1:total_time/2,:)/2; end
% alter source conditions
if 0; source(total_time/2:end,:)=source(total_time/2:end,:)/2;end
% set fixed wind climate

% Ct (=concentration in transport)
Ct = zeros(total_time/dt,L_dom/dx+1);
% Ca (=concentration available on bed)
Ca = zeros(total_time/dt,L_dom/dx+1);
Ca(1,:)=source(1,:);

% Cu =(concentration based on wind speed)
Cu = zeros(total_time/dt,L_dom/dx+1);

[Ct,Ca,Cu,Ccap_index] = model_core(u_w,U_th,source,dx,dt,total_time,T,VS,z);

% plot_results(Ct,Ca,u_w,source,dt,1)

%% plotting proves to be a bit difficult.
% Lets try the panel package to improve the plotting.
addpath panel
figure(5)
p = panel();

% let's start with three columns, both halve
p.pack('h', [1/2 1/2])

% then let's pack 3 rows into the first column
p(2).pack([1/6 1/3 1/3 1/6]);

% then let's pack 2x6 rows into the seceond column
p(1).pack('h', [1/2 1/2]);

p(1,1).pack(6,1);
p(1,2).pack(6,1);

p.de.margin = 2;

p.margintop = 10;

p(1,2).marginleft = 15;
p(2).marginleft = 25;
p(2).margintop = 25;

p(2,1).marginbottom = 10;
p(2,2).marginbottom = 10;
p(2,3).marginbottom = 10;

times = [15 30 45 60 75 90];

p(2,2).select()
plot(dt:dt:total_time,u_w(:,101),'linewidth',2,'Color',[0.4 0.4 0.4])
ylim([0 11])
box on
vline(times,'k:')
ylabel('Wind [m/s]')
set(gca, 'xtick', [])

p(2,3).select()
plot(dt:dt:total_time,u_w(:,101).*Ct(:,101),'Color',[0.4 0.4 0.4],'linewidth',2)
hold all
plot(dt:dt:total_time,1.5e-4*(u_w(:,101)-U_th).^3,'k--')%,'Color',[0.4 0.4 0.4]
ylim([0 0.02])
vline(times,'k:')
box on
ylabel('Q at x=101m [kg/ms]')
xlabel('Time [s]')



set(gca, 'xtick', times)


for i = 1:6
    p(1,1,i,1).select()
    plot(Ct(times(i)/dt,:).*u_w(times(i)/dt,:),'linewidth',2,'Color',[0.4 0.4 0.4])
%     hline(u_w(times(i),1)^3/10^3,'--k')
    hline(1.5e-4*((u_w(times(i)/dt,1)-U_th).^3),'--k')
    
    ylim([0 0.03])
    xlim([0 110])
    box on
    text(30,0.026,['t = ' num2str(times(i))])
    
    if i~=6
        set(gca, 'xtick', [])
    end
    
    p(1,2,i,1).select()
    plot(Ca(times(i)/dt,:),'linewidth',2,'Color',[0.4 0.4 0.4])
    ylim([0 1.1*max(max(Ca))])
    xlim([0 110])
    box on
    text(30,0.95*max(max(Ca)),['t = ' num2str(times(i))])
    set(gca,'ytick',[0 0.01])
    if i~=6
        set(gca, 'xtick', [])
    end
    
end

p(1).xlabel('Distance - x [m]');
p(1,1).ylabel('Q [kg/s/m]')
p(1,2).ylabel('Erodible sediment at the bed (Se) [Kg/m^2]')

p(1).title('Spatial series')
p(2,2).title('Time series')


