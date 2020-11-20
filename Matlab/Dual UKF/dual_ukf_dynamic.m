 %% Estimation of Tire Co-efficients and Mass using Dynamic model
close all

clear all
%loading control input data
load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

dt= 0.1;
pre_computed_params = [1.56,1.385];
orig_params = [1.5,1.5];
%dynamic states = vx vy yaw_rate 
%observing vx yaw_rate
x = [1;0;0.2];
%tire stiffness co-efficients
param = [1700; 15e4; 5e4];

% The system model for both UKF's
xsbar = [0.5;0;0.1]; % initial belief of state
xpbar = [1500; 18e4; 4e4]; % initial belief of parameter

Ps = diag([1, 1, 0.2]); %initial belief uncertainity
Pp = diag([200 2e4 1.5e4]);
Rs = diag([1.2^2 1.2^2 0.45^2]); %process noise
Rp = diag([5^2 5.75e2^2 2.75e2^2]);

Q = 10*diag([0.02^2 0.02^2 0.025^2]); %measurement noise
Hs = eye(3);
gamma = 0.998;
n = size(u1,1);
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;
X(:,1) = x;
Ppp=[Pp(1,1);Pp(2,2)];
err_s = [];
err_p = [];

disp('Estimating parameters from dynamic model')
time = [];
current_time = 0;

for loop = 1:n
    
    time = [time;current_time];
    current_time = current_time+dt;
    
    u = u1(loop,:);
    x = dynamic_step(x,u,param,orig_params,dt,0.01);
    z = Hs*x + 0.02*randn(3,1);
    %% Prediction
    %%Parameter
    %UKF defining tuning variables and weights
    L_p=numel(xpbar);
    alpha_p=1e-2; %tune here.. for parameters
    ki_p=2;
    beta_p=2;
    lambda_p=alpha_p^2*(L_p+ki_p)-L_p;
    c_p=L_p+lambda_p;
    Wm_p=[lambda_p/c_p 0.5/c_p+zeros(1,2*L_p)];           %weights for means
    Wc_p=Wm_p;
    Wc_p(1)=Wc_p(1)+(1-alpha_p^2+beta_p);               %weights for covariance
    c_p=sqrt(c_p);
    
    %UKF defining sigma points Xt-1, Step 2 from book
    Pp = (Pp + Pp')/2;
    Xp=sigmas(xpbar,Pp,c_p);                            %sigma points around x
    
    %parameter prediction using unscented transform
    xphat = zeros(L_p,1);
    for l=1:size(Xp,2)
        Yp(:,l)=Xp(:,l); %The map of parameter is an identity
        xphat=xphat+Wm_p(l)*Yp(:,l);
    end
    Y1p=Yp-xphat(:,ones(1,size(Xp,2)));
    Pphat=Y1p*diag(Wc_p)*Y1p'+Rp;
    
    %%Taking sigma points around the predicted state %State 6: from book
    Pphat = (Pphat + Pphat')/2;
    Xphat=sigmas(xphat,Pphat,c_p);


    %%State
    %UKF defining tuning variables and weights
    L_s=numel(xsbar);
    alpha_s=1e-3;
    ki_s=0;
    beta_s=2;
    lambda_s=alpha_s^2*(L_s+ki_s)-L_s;
    c_s=L_s+lambda_s;
    Wm_s=[lambda_s/c_s 0.5/c_s+zeros(1,2*L_s)];           %weights for means
    Wc_s=Wm_s;
    Wc_s(1)=Wc_s(1)+(1-alpha_s^2+beta_s);               %weights for covariance
    c_s=sqrt(c_s);
    
    %UKF defining sigma points Xt-1
    Ps = (Ps + Ps')/2;
    Xs=sigmas(xsbar,Ps,c_s);                            %sigma points around x
    xsbar_prev = xsbar;
    %state prediction using unscented transform
    xshat = zeros(L_s,1);
    for l=1:size(Xs,2)%check if it is 2N or not
        Ys(:,l)=dynamic_step(Xs(:,l),u,xphat,pre_computed_params,dt,0); %The map of state
        xshat=xshat+Wm_s(l)*Ys(:,l);
    end
    Y1s=Ys-xshat(:,ones(1,size(Xs,2)));
    Pshat=Y1s*diag(Wc_s)*Y1s'+Rs;
    
    %%Taking sigma points around the predicted state %State 6: from book
    Pshat = (Pshat + Pshat')/2;
    Xshat=sigmas(xshat,Pshat,c_s);

    
    %% Correction
    %%UKF on measurement for state
    
    %measurement prediction using unscented transform
    %parameter prediction using unscented transform, Step 3,4,5: from book
    zshat = zeros(L_s,1);
    for l=1:size(Xs,2)%check if it is 2N or not
        Yz(:,l)=Xshat(:,l); %The map of state
        zshat=zshat+Wm_s(l)*Yz(:,l);
    end
    
    Y1z=Yz-zshat(:,ones(1,size(Yz,2)));
    Pzhat=Y1z*diag(Wc_s)*Y1z'+Q;
    %%Evaluating cross covariance. ilne 10 from book
    Psz=Y1s*diag(Wc_s)*Y1z';
       
    %%State correction
     
    Ks=Psz*inv(Pzhat);
    xsbar=xshat+Ks*(z-zshat);                              %state update
    Ps=Pshat-Ks*Pzhat*Ks';                                %covariance update

    
    %Hope it is correct.. !!
    
    %%UKF on measurement for parameters
    zphat = zeros(L_s,1);
    for l=1:size(Xp,2)%check if it is 2N or not
        Yzp(:,l) = dynamic_step(xsbar_prev,u,Xp(:,l),pre_computed_params,dt,0); %Notice that I am generating expected measurement for each sigma point of parameter distribution
        zphat=zphat+Wm_p(l)*Yzp(:,l);
    end
    Y1zp=Yzp-zphat(:,ones(1,size(Xp,2)));
    Pzhat_2=Y1zp*diag(Wc_p)*Y1zp'+Q;
    Ppz = Y1p*diag(Wc_p)*Y1zp';  %critical point, generating cross covariance
    
    Kp=Ppz*inv(Pzhat_2);
    xpbar=xphat+Kp*(z-zphat);                              %state update
    Pp=Pphat-Kp*Pzhat_2*Kp';                            %covariance update

    %disp(Kp*Ppz')
    Rp = Rp*gamma;
    Xsbar(:,loop+1) = xsbar;
    Xpbar(:,loop+1) = xpbar;
    X(:,loop+1) = x;
    Ppp(:,loop+1) = [Pp(1,1);Pp(2,2)];
    s_err = xsbar - x;
    err_s = [err_s s_err];
    p_err = xpbar - [1700;15e4;5e4];
    err_p = [err_p p_err];
end
figure,
subplot(3,1,1)
plot(Xsbar(1,:),'r')
hold on
plot(X(1,:),'b')
subplot(3,1,2)
plot(Xsbar(2,:),'r')
hold on;
plot(X(2,:),'b')
subplot(3,1,3)
plot(Xsbar(3,:),'r')
hold on;
plot(X(3,:),'b')

figure
subplot(3,1,1)
plot(err_s(1,:))
subplot(3,1,2)
plot(err_s(2,:))
subplot(3,1,3)
plot(err_s(3,:))

figure,
subplot(3,1,1)
plot(Xpbar(1,:),'k')
subplot(3,1,2)
plot(Xpbar(2,:),'m')
subplot(3,1,3)
plot(Xpbar(3,:),'c')

% figure
% subplot(3,1,1)
% plot(err_p(1,:))
% subplot(3,1,2)
% plot(err_p(2,:))
% subplot(3,1,3)
% plot(err_p(3,:))

disp('Done')
Xpbar(:,loop+1)
mean(abs(err_s),2)
% mean(abs(err_p(:,500,end)),2)
% 
% figure
% plot(Ppp(1,:), 'r')
% hold on
% plot(Ppp(2,:),'b')