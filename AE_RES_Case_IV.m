% 3rd Power concentration

%% now lets assume some variability in wind.
% therefore we have dc(t)/dt+dc(x)u(x)/dx=source(x)
% now we need a time varying concentration. Let's assume for the moment
% that the spatial wind profile migrates with constant speed in the
% direction of the wind and that the concentration migrates likewhise.


clear;close all
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
    %     f_mean = 10; %m/s
    f_sigma = 2; %m/s
    %     f_sigma = 3; %m/s
    
    
    % length
    l_mean = 4; %s
    l_sigma = 2; %s
    
    u_w = GenWind(f_mean,f_sigma,l_mean,l_sigma,total_time,dt);
    u_w = repmat(u_w,1,L_dom/dx+1);
    
else
    % Try if we can reuse the previously correct dataset
    disp('loading wind data from disc')
    %     load wind_sep.mat
    %     load sep_9okt.mat
    load sep_10oktc.mat
    
end



% source magnitude(s)
s = [0.1:0.1:5 6:2:16]*1e-4;

% Threshold velocity
U_th = 4; %[m/s]

% adaptation timescale
T=0.5; %[s]
VS=1;
z=[];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% initiate empty matrices
% Ct (=concentration in transport)
Ct = zeros(total_time/dt,L_dom/dx+1);
% Ca (=concentration available on bed)
Ca = zeros(total_time/dt,L_dom/dx+1);
% Cu =(concentration based on wind speed)
Cu = zeros(total_time/dt,L_dom/dx+1);
% source
source=zeros(length(s),total_time/dt,L_dom/dx+1);

for k=1:length(s)
    source(k,:,2:L_dom*alpha/dx)= s(k)*dt*dx;
end
%
% CT101
Ct101= zeros(total_time/dt,length(s));




R_lin = zeros(1,length(s));
R_3rd = zeros(1,length(s));


for k = 1:length(s)

    [Ct,Ca,Cu,Ccap_index] = model_core(u_w,U_th,squeeze(source(k,:,:)),dx,dt,total_time,T,VS,z);
    
   
    
    % some resampling to reduce matrices in a similar way data is often
    % treathed.
            samples = 1; % per 5 sec
        u_1hz = u_w(1:samples/dt:end,end);
    q_1hz = (VS*u_w(1:samples/dt:end,end).*Ct(1:samples/dt:end,end));
    
        % calculate averages over a window (this is also a similar as common data processing. 
        window = 10 ; %seconds
    for i = 1:floor(length(u_1hz)/window)
        u_fit(k,i) = mean(u_1hz(((i-1)*window)+1:i*window));
        q_fit(k,i) = mean(q_1hz(((i-1)*window)+1:i*window));
    end
    
%     Here we fit our linear and 3rd power functions
    [V_th3(k),cub_coef(k),R_3rd(k)] = powerfit(u_fit(k,:),q_fit(k,:),0,'3rd');
    [V_th1(k), lin_coef(k),R_lin(k)] = powerfit(u_fit(k,:),q_fit(k,:),0,'lin');

    
    %% now we define Q at the end of the domain for this source magnitude
    Ct101(:,k) = Ct(:,end);
    
end

%% now we plot Figure 8 using the panel package
figure(121);close;figure(121)
p = panel();

% lets define 2 rows
p.pack('v', [1/100 (2/5-1/100) 1/5 1/5 1/5]);

% in the first row there are three plots
p(2).pack('h', [1/3 1/3 1/3])
p(3).pack('h', [1/10 8/10 1/10])
p(4).pack('h', [1/10 8/10 1/10])
p(5).pack('h', [1/10 8/10 1/10])


p.de.margin = 2;
p(2).margintop =2;
p(2,2).marginleft =10;
p(2,3).marginleft =10;


p(3).margintop =15;
p(3,2).marginleft =10;
p(3,2).marginright =10;

p(4).margintop =5;
p(4,2).marginleft =10;
p(4,2).marginright =10;

p(5).margintop =5;
p(5,2).marginleft =10;
p(5,2).marginright =10;

out = [1.5 5 10]*1e-4;

for i=1:3;
    index = find(abs(s-out(i))<0.0000001);
    
    figure(121)
    p(2,i).select()
    hold all
    
    plot([0 V_th1(index):0.01:max(u_w(:,end))],...
        [0 lin_coef(index)*([V_th1(index):0.01:max(u_w(:,end))]-V_th1(index))],...
        'k-.','linewidth',2,'color',[0.5 0.5 0.5])
    
    plot([0 V_th3(index):0.01:max(u_w(:,end))],...
        [0 cub_coef(index)*([V_th3(index):0.01:max(u_w(:,end))]-V_th3(index)).^3],...
        'k--','linewidth',2,'color',[0 0 0])
    
    legend(['R^2 = ' num2str(round(R_lin(index)*100)/100)],['R^2 = ' num2str(round(R_3rd(index)*100)/100)],'Location','NW')
    
    
    
    plot(u_fit(index,:),q_fit(index,:),'.','color',[0.8 0.8 0.8])
    
    plot([0 U_th:0.01:max(u_w(:,end))],...
        [0 1.5e-4*([U_th:0.01:max(u_w(:,end))]-U_th).^3],...
        'k-','linewidth',1)
    
    plot([0 V_th3(index):0.01:max(u_w(:,end))],...
        [0 cub_coef(index)*([V_th3(index):0.01:max(u_w(:,end))]-V_th3(index)).^3],...
        'k--','linewidth',2,'color',[0 0 0])
    
    
    plot([0 V_th1(index):0.01:max(u_w(:,end))],...
        [0 lin_coef(index)*([V_th1(index):0.01:max(u_w(:,end))]-V_th1(index))],...
        'k-.','linewidth',2,'color',[0.5 0.5 0.5])
    
    legend(['R^2 = ' num2str(round((R_lin(index))*100)/100)],['R^2 = ' num2str(round((R_3rd(index))*100)/100)],'Location','NW')
    %
    if i==2
        title({'B. Supply \propto wind capacity'})
        xlabel('Wind speed [m/s]')
        
    elseif i==3
        title('C. Supply >> wind capacity')
    elseif i==1
        title('A. Supply << wind capacity')
        ylabel('Q [Kg/s/m]')    
    end
    ylim([0 0.08])
    xlim([0 12])
    
    box on

end
%
p(3,2).select()
plot(s*L_dom*alpha*total_time/1e2,(R_3rd),'k-.','linewidth',2)
hold all
plot(s*L_dom*alpha*total_time/1e2,(R_lin),'k-','linewidth',2)
ylim([0 1])

vline(out*L_dom*alpha*total_time/1e2,'k-')

set(gca, 'xtick', []);
ylabel('R^2 [-]');
legend('Cubic fit','Linear fit','Location','SE')
box on

text(out(1)*L_dom*alpha*total_time/1e2+0.01,0.1 ,'A')
text(out(2)*L_dom*alpha*total_time/1e2+0.01,0.1,'B')
text(out(3)*L_dom*alpha*total_time/1e2+0.01,0.1,'C')

p(4,2).select()
plot(s*L_dom*alpha*total_time/1e2,lin_coef,'k','linewidth',2)
hold all
plot(s*L_dom*alpha*total_time/1e2,cub_coef*1e2,'k-.','linewidth',2)
set(gca, 'xtick', []);
ylabel([{'A_{cub} [10^2 kg s^2/m^4]'; 'A_{lin} [kg/m^2]'}]);%;{'(from linear fit)'}])
box on
vline(out*L_dom*alpha*total_time/1e2,'k-');

text(out(1)*L_dom*alpha*total_time/1e2+0.01,1.8 ,'A')
text(out(2)*L_dom*alpha*total_time/1e2+0.01,1.8,'B')
text(out(3)*L_dom*alpha*total_time/1e2+0.01,1.8,'C')

p(5,2).select()
plot(s*L_dom*alpha*total_time/1e2,V_th1,'k','linewidth',2)
hold all
plot(s*L_dom*alpha*total_time/1e2,V_th3,'k-.','linewidth',2)
ylim([0 8])

ll=hline(U_th,':k');
xlabel('\int \int S_s dt dx [10^2 kg/m]')
ylabel([{'V_{th_{cub}} [m/s]';'V_{th_{lin}} [m/s]'} ]);%;{'(from linear fit)'}])
vline(out*L_dom*alpha*total_time/1e2,'k-')

box on

if 0 % save figure
p.export('Figure8.eps', '-w150', '-h200','-oeps');
end


%% Here we plot Figure 9
figure
plot(s*L_dom*alpha*total_time/1e2,(mean(u_w(:,end))-V_th1).*lin_coef*3600/1e2,'k','linewidth',2,'DisplayName','$Q=A_{lin}(\overline{u_w}-u_t)$')
hold all

sup_Est1(k) = sum((u_w(((u_w(:,end)-V_th1(k))>0),end)-V_th1(k))*lin_coef(k)*dt);

plot(s*L_dom*alpha*total_time/1e2,s*L_dom*alpha*total_time/1e2,'k--','linewidth',1,'DisplayName','$\int Q dt =\int\int S_s dt dx$')
xlabel('\int \int S_s dt dx [10^2 kg/m]')
ylabel([{'\int Q dt [10^2 kg/m]'}]);%;{'(from linear fit)'}])
ylim([0 0.8])
xlim([0 0.8])
box on

plot(s*L_dom*alpha*total_time/1e2,sum(Ct101(:,:).*repmat(u_w(:,end),1,length(s)))*dt/1e2,'k-','DisplayName','Modelled\ $Q$')
plot([0 100],[sum((1.5e-4*(u_w(:,end)-U_th).^3))*dt/1e2 sum((1.5e-4*(u_w(:,end)-U_th).^3))*dt/1e2],':k','linewidth',1,'DisplayName','$Q_u$')
l= legend('show','location','EO');
set(l,'Interpreter','latex','Fontsize',8,'Fontweight','Bold');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.634517 6.34517 12 6])
