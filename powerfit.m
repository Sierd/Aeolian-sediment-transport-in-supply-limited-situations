function [V_th,C,R_sq] = powerfit(u,q,plots,mode)
% This function is made in order to generate a generic fitting tool for
% third power functions

%{
to test the following script can be run in command prompt.
    close all
    clear
    u = 1:0.01:10;
    q = 2*(u-5).^3;
    q2 = 2*(u-5);
    q(q<0) = 0;
    q2(q2<0) = 0;

%     q(q>0) = q(q>0)+500;
    [V_th,c,P_corr] = powerfit(u,q,1,'lin')
    [V_th,c,P_corr] = powerfit(u,q,1,'3rd')
    
    [V_th,c,P_corr] = powerfit(u,q2,1,'lin')
    [V_th,c,P_corr] = powerfit(u,q2,1,'3rd')

%}

if ~exist('plots');
    plots = 0;
end

if mode=='3rd'
    %% let's iterate threshold velocities.
    vt = 0:0.1:max(u);
    SS_tot = sum((q-mean(q)).^2);
 
    c = 0:0.1e-6:1e-4;
    R_sq_all=ones(length(c),length(vt));
    
    for tel_VT =1:length(vt)
        %         for tel_C = 1:length(c)
        mod_data = c'*(u-vt(tel_VT)).^3;
        mod_data(mod_data<0)=0;
        % SS_reg = sum((polyval(p_line,x)-mean(y)).^2);
        SS_err = sum((repmat(q,length(c),1)-mod_data).^2,2);
        R_sq_all(:,tel_VT) = 1-SS_err/SS_tot;
        %         end
    end
    [R_sq,V_loc] = max(max(R_sq_all,[],1));
    [R_sq,C_loc] = max(max(R_sq_all,[],2));
    
    V_th = vt(V_loc);
    
    %         C_prev=C;
    C =  c(C_loc);
    
    
    
    if plots
        %         figure;
        % here we plot the data
        plot(u,q,'.','color',[0.8 0.8 0.8])
        hold all
        % here we plot the fitted model
        mod_data = C*((1:0.1:max(u))-V_th).^3;
        mod_data(((1:0.1:max(u))-V_th)<0)=0;
        plot(1:0.1:max(u),mod_data,'--','color',[0 0 0],'linewidth',2)
    end
    
elseif mode=='lin'
    vt = 0:0.1:max(u);
    %     c = 0:0.001:1;
    %     c = 0:0.1e-4:5e-3;
    
    SS_tot = sum((q-mean(q)).^2);
    
    
    
    c = 0:0.1e-4:1e-2;
    R_sq_all=ones(length(c),length(vt));
    
    
    
    for tel_VT =1:length(vt)
        %         for tel_C = 1:length(c)
        mod_data = c'*(u-vt(tel_VT)).^1;
        mod_data(mod_data<0)=0;
        % SS_reg = sum((polyval(p_line,x)-mean(y)).^2);
        SS_err = sum((repmat(q,length(c),1)-mod_data).^2,2);
        R_sq_all(:,tel_VT) = 1-SS_err/SS_tot;
        %         end
    end
    [R_sq,V_loc] = max(max(R_sq_all,[],1));
    [R_sq,C_loc] = max(max(R_sq_all,[],2));
    V_th = vt(V_loc);
    C =  c(C_loc);
    
    
    
    
    
    
    if plots
        %         figure;
        % here we plot the data
        plot(u,q,'*','color',[0 0 0])
        hold all
        % here we plot the fitted model
        mod_data = C*((1:0.1:max(u))-V_th).^1;
        mod_data(((1:0.1:max(u))-V_th)<0)=0;
        plot(1:0.1:max(u),mod_data,'-.','color',[0.5 0.5 0.5],'linewidth',2)
    end
    
    
    
    
end



