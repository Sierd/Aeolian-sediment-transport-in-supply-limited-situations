% AE_RES Test Case II
%%%%%%%%%%%%%%%%%%%%%%%%%%%   INPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation time
total_time = 3600; %[s]


% length of domain
L_dom = 100;% [m]
% relativel length of supply zone
alpha = 0.2;



% numerical stuff
dx=1;dt=0.05;
DT=dt;

% generate random varying wind
if 0
    disp('creating new wind data')
    % force
    f_mean = 7; %m/s
    f_sigma = 2.5; %m/s
    
    % length
    l_mean = 4; %s
    l_sigma = 4; %s
    
    u_w = GenWind(f_mean,f_sigma,l_mean,l_sigma,total_time,dt);
    u_w = repmat(u_w,1,L_dom/dx+1);
    
else
    % Try if we can reuse the previously correct dataset
    disp('loading wind data from disc')
    %     load wind_sep.mat
    load wind.mat
    %     load case3sep2.mat
end


% source magnitude(s)
src = [1.5e-4]; %[Kg/m2s]
% src = [1.4e-3]; %[Kg/m2s] lekker veel !!!
% src = [10]; %[Kg/m2s]
% Threshold velocity
U_th = 4; %[m/s]



% adaptation timescale
T=0.5; %[s]



%%%%%%%%%%%%%%%%%%%% EXTRA PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%
VS=1;% [m ???]
z=0.1;% [m]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
source=zeros(total_time/dt,L_dom/dx+1);
source(:,2:L_dom*alpha/dx)= src*dt*dx;

x_line = 0:dx:L_dom;
t_line = 0:dt:total_time-dt;

% Ct (=concentration in transport)
Ct = zeros(total_time/dt,L_dom/dx+1);
% Ca (=concentration available on bed)
Ca = zeros(total_time/dt,L_dom/dx+1);
Ca(1,:)=source(1,:);

% Cu =(concentration based on wind speed)
Cu = zeros(total_time/dt,L_dom/dx+1);

tic
[Ct,Ca,Cu,Ccap_index] = model_core(u_w,U_th,source,dx,dt,total_time,T,VS,z);
toc

%  plot_results_case3(Ct,Ca,u_w,source,dt,1)



%% FITTING !!
% Lets plot q and v averaged over 1 second.
samples = 1; % second interval
% samples = dt; % all

u_1hz = u_w(1:samples/dt:end,end);
q_1hz = (VS*u_w(1:samples/dt:end,end).*Ct(1:samples/dt:end,end));

clear u_fit q_fit

window = 10 ; %seconds
% window = 1 ; %all
for i = 1:floor(length(u_1hz)/window)
    u_fit(i) = mean(u_1hz(((i-1)*window)+1:i*window));
    q_fit(i) = mean(q_1hz(((i-1)*window)+1:i*window));
end


figure(7)
[V_th3,c_3,P_3] = powerfit(u_fit,q_fit,1,'3rd');

plot([0 U_th:0.1:15],1.5e-4*(([U_th U_th:0.1:15])-U_th).^3,'k','linewidth',2)

init = (1.5e-4*(u_fit-U_th).^3);
init(init<0)=0;

[rsq_ini] =1- sum((q_fit-init).^2)/sum((q_fit-mean(q_fit)).^2);


[V_th1,c_1,P_1] = powerfit(u_fit,q_fit,1,'lin');
% [V_th1NEW,c_1NEW,P_1NEW] = powerfit_new(u_fit,q_fit,1,'lin')


ylim([0 max(q_fit)])


xlabel('Wind speed [m/s]')
ylabel('Q [kg/s] at x=101m')
% ylim([0 1.1])
%  xlim([0 10])
% xlim([0 12])
set(gcf,'paperunits','centimeters','PaperPosition',[0.634517 6.34517 7 5])



%%
addpath panel


figure(6)
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

times = [600:600:3600];


p(2,2).select()
plot(u_w(1:1/dt:end,101),'linewidth',2,'Color',[0.4 0.4 0.4])
hline(mean(mean(u_w)),'--k')
ylim([0 15])
xlim([0 3650])

box on
vline(times,'k:')
ylabel('Wind [m/s]')
set(gca, 'xtick', [])

p(2,3).select()
%     plot(u_w(1:1/dt:end,101).*Ct(1:1/dt:end,101)/max(max(u_w))^3,'Color',[0.4 0.4 0.4],'linewidth',2)
plot(VS*u_w(1:1/dt:end,101).*Ct(1:1/dt:end,101),'Color',[0.4 0.4 0.4],'linewidth',2)
% hold all
% plot(u_w(:,101).^3/max(max(u_w))^3,'k--')%,'Color',[0.4 0.4 0.4]
xlim([0 3650])
vline(times,'k:')
box on
ylabel('Q at x=101m [kg/ms]')
xlabel('Time [s]')


set(gca, 'xtick', times)


for i = 1:6
    p(1,1,i,1).select()
    %     plot(Ct(times(i),:).*u_w(times(i),:)/max(max(u_w-U_th))^3,'linewidth',2,'Color',[0.4 0.4 0.4])
    plot(VS*Ct(times(i),:).*u_w(times(i),:),'linewidth',2,'Color',[0.4 0.4 0.4])
    %         hline(u_w(times(i),1)^3/max(max(u_w))^3,'--k')
    if 1.5e-4*(u_w(times(i),1)-U_th)^3>0
        hline(1.5e-4*(u_w(times(i),1)-U_th)^3,'--k')
    else
        hline(0,'--k')
    end
    ylim([0 0.025])
    %         ylim(1.1*ylim)
    xlim([0 110])
    box on
    text(30,0.85*max(ylim),['t = ' num2str(times(i))])
    if i~=6
        set(gca, 'xtick', [])
    end
    
    p(1,2,i,1).select()
    plot(Ca(times(i),:)*1e3,'linewidth',2,'Color',[0.4 0.4 0.4])
    ylim([0 1.1*max(max(Ca(times,:)*1e3))])
    xlim([0 110])
    box on
    text(30,0.85*max(ylim),['t = ' num2str(times(i))])
    if i~=6
        set(gca, 'xtick', [])
    end
    
end

p(1).xlabel('Distance - x [m]');
p(1,1).ylabel('Q [kg/s/m]')
p(1,2).ylabel('Erodible sediment at the bed (Se) [10^3 Kg/m^2]')

p(1).title('Spatial series')
p(2,2).title('Time series')
% p.refresh();
if 0
    p.export('Case_3sep.eps', '-w170', '-h150','-oeps');
end


