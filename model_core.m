% This function is the core of the 1st order upwind model.
% Ct is the Sediment concentration in transport
% u_w is the wind timeseries
% Ca is the sediment concetration at the bed
% Cu is the sediment concentration as a function of wind speed
% Cth is the threshold sediment concentration
% source is the source area
% dx and dt are descretization constants.
% T is a pick up timescale

function [Ct,Ca,Cu,Ccap_index] = model_core(u_w,U_th,source,dx,dt,~,T,VS,~)


% test for stability
if dx/dt<max(u_w)
    error('dx/dt < max(u_w), please adapt numerical parameters')
end

% Ct (=concentration in transport)
Ct = zeros(size(source));
Ct1 = zeros(size(source));
Ct2 = zeros(size(source));

% Ca (=concentration available on bed)
Ca = zeros(size(source));
Ca(1,:)=source(1,:);

% Cu =(concentration based on wind speed)
Cu = zeros(size(source));

Ccap_index = zeros(size(source));

%first calculate sediment transport rate according to bagnold. This can
%take place outside the time loop which is better for speed
Cu(:,:) = 1.5e-4*((u_w(:,:)-U_th).^3)./(u_w(:,:)*VS);
%         Cu(t,:) = (u_w(t,i)-U_th)^2;
%                 Cu(t,:) = (u_w(t,i)-U_th)*u_w(t,i);
%                 Cu(t,:) = ((u_w(t,i)-U_th)*(u_w(t,i)+U_th)^2)/u_w(t,i);
%remove all values smaller than 0
Cu(Cu<0)=0;


% h = waitbar(0,'Model runs');

for t = 2:size(Ct(:,1),1)
    %     waitbar(t/size(Ct(:,1),1))
    %upwind boundary 1 cell
    %     i=1;
    
    % calculate Cu (=wind driven sediment concentration)
    
    % impose Ca (=concentration available on bed) in source area
%     Ca(t,:)=Ca(t-1,:)+source(t,:)/dx;
    
    %%% lets try a matrix solver
    % we have u_w, Ct and Ca
    % we assume Ct(1) = 0
    %           Ca(1) = 0
    
    % This first solution is valid when enough sediment is available at the bed
    Ct1(t,2:end) = ((-VS*u_w(t-1,1:end-1).*(Ct(t-1,2:end)-Ct(t-1,1:end-1))/dx)*dt...
        +Ct(t-1,2:end)+Cu(t-1,2:end)/(T/dt))/(1+1/(T/dt));
    
    % create index for 2nd solution
    index =  (Cu(t-1,:)-Ct1(t,:))/(T/dt) > Ca(t-1,:);
    % This second solution is valid when sediment availability at the bed is not enough
    Ct2(t,2:end) = (-VS*u_w(t-1,1:end-1).*(Ct(t-1,2:end)-Ct(t-1,1:end-1))/dx)*dt...
        +Ct(t-1,2:end)+Ca(t-1,2:end)/(T/dt);
    
    Ct(t,~index)=Ct1(t,~index);
    Ct(t,index)=Ct2(t,index);
   
    Ca(t,index)=Ca(t-1,index)+source(t,index)/dx-Ca(t-1,index)/(T/dt);
    Ca(t,~index)=Ca(t-1,~index)+source(t,~index)/dx-(Cu(t-1,~index)-Ct(t,~index))/(T/dt);
%     Ca(t,~index)=Ca(t,~index)-(Cu(t-1,~index)-Ct1(t-1,~index))/(T/dt);
    
    %      Ccap_index(t,i) = 1;
 end
% close(h)

