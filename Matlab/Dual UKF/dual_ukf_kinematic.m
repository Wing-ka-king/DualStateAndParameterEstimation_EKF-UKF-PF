%% Dual UKF based state and parameter estimation for autonomous vehicle control
% 
clc
clear all
close all

%% Estimation of COM using Kinematic model - lf, lr
dt= 0.1;

%params = [lf,lr]
%states = [x,y,phi,v];
%u = [acc,steering_angle]

x = [0;0;0;0];
params = [1.5; 1.5];

%The system model for both UKF's
xsbar = [0; 0; 0; 0]; % initial belief
xpbar = [4; 4]; % initial belief

Ps = diag([1, 1, 0.2, 2]); %initial belief uncertainity
Pp = diag([3, 3]);
Rs = diag([0.5^2 0.5^2 0.2^2 1^2]); %process noise
Rp = diag([0.25^2 0.25^2]);
Q = diag([1^2 1^2 0.1^2 0.1^2]); %measurement noise

n = 1000; %number of iteration

X(:,1) = x;
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;
sigmaP = Pp(:);
sigmaS = Ps(:);
errpose = [];

sim_noise = 0.1;
gamma = 0.97;
disp('Estimating parameters from kinematic model')
for k = 1:n
    
    u = [2.5*cos(0.01*k),0.15*sin(0.01*k)];
    %u = [7.5e5,0.1];
    x = kinematic_step(x,params,u,dt,sim_noise);
    z = x + 0.1*randn(4,1); %observing all state in simulation mode
    
    %% Prediction
    %% Parameter
    %UKF defining tuning variables and weights
    L_p=numel(xpbar);
    alpha_p=1e-1; %tune here.. for parameters
    
    ki_p=2;
    beta_p=2;
    lambda_p=alpha_p^2*(L_p+ki_p)-L_p;
    c_p=L_p+lambda_p;
    Wm_p=[lambda_p/c_p 0.5/c_p+zeros(1,2*L_p)];           %weights for means
    Wc_p=Wm_p;
    Wc_p(1)=Wc_p(1)+(1-alpha_p^2+beta_p);               %weights for covariance
    c_p=sqrt(c_p);
    
    %UKF defining sigma points Xt-1  %% Step 2: from book
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
    
    %% State
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
    
    %UKF defining sigma points Xt-1 %% Step 2: from book
    Ps = (Ps + Ps')/2;
    Xs=sigmas(xsbar,Ps,c_s);                            %sigma points around x
    xsbar_prev = xsbar;
    %parameter prediction using unscented transform, Step 3,4,5: from book
    xshat = zeros(L_s,1);
    for l=1:size(Xs,2)%check if it is 2N or not
        Ys(:,l)=kinematic_step(Xs(:,l),xphat,u,dt,0); %The map of state
        xshat=xshat+Wm_s(l)*Ys(:,l);
    end
    Y1s=Ys-xshat(:,ones(1,size(Xs,2)));
    Pshat=Y1s*diag(Wc_s)*Y1s'+Rs;
    
    %%Taking sigma points around the predicted state %State 6: from book
    Pshat = (Pshat + Pshat')/2;
    Xshat=sigmas(xshat,Pshat,c_s);
    
    %% Correction
    %% State
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
    
    
    %%UKF on measurement for parameters
    zphat = zeros(L_s,1);
    for l=1:size(Xphat,2)%check if it is 2N or not
        Yzp(:,l) = kinematic_step(xsbar_prev,Xphat(:,l),u,dt,0); %Notice that I am generating expected measurement for each sigma point of parameter distribution
        zphat=zphat+Wm_p(l)*Yzp(:,l);
    end
    Y1zp=Yzp-zphat(:,ones(1,size(Xp,2)));
    Pzhat_2=Y1zp*diag(Wc_p)*Y1zp'+Q;
    Ppz = Y1p*diag(Wc_p)*Y1zp';  %critical point, generating cross covariance
    
    Kp=Ppz*inv(Pzhat_2);
    xpbar=xphat+Kp*(z-zphat);                              %state update
    Pp=Pphat-Kp*Pzhat_2*Kp';                                %covariance update
    
    Xsbar(:,k+1) = xsbar;
    Xpbar(:,k+1) = xpbar;
    X(:,k+1) = x;
    Rp = gamma*Rp;
    
    sigmaP = [sigmaP Pp(:)];
    sigmaS = [sigmaS Ps(:)];
    err = xpbar - [1.5;1.5];
    errpose = [errpose err];
end
% figure
% subplot(4,1,1)
% plot(Xsbar(1,:),'b')
% hold on;
% plot(X(1,:),'r')
% subplot(4,1,2)
% plot(Xsbar(2,:),'b')
% hold on;
% plot(X(2,:),'r')
% subplot(4,1,3)
% plot(Xsbar(3,:),'b')
% hold on;
% plot(X(3,:),'r')
% subplot(4,1,4)
% plot(Xsbar(4,:),'r')
% hold on;
% plot(X(4,:),'b')
% 
% 
% figure
% subplot(4,1,1)
% plot(sigmaS(1,:))
% subplot(4,1,2)
% plot(sigmaS(2,:))
% subplot(4,1,3)
% plot(sigmaS(3,:))
% subplot(4,1,4)
% plot(sigmaS(4,:))



[Xpbar(1,n), Xpbar(2,n)],
diff = [abs(Xpbar(1,n)-1.5)/1.5*100, abs(Xpbar(2,n)-1.5)/1.5*100]

if (diff(1) && diff(2) < 5)
    figure,
    subplot(2,1,1)
    plot(Xpbar(1,:),'r')
    subplot(2,1,2)
    plot(Xpbar(2,:),'b')
    
    figure
    subplot(2,1,1)
    plot(sigmaP(1,:))
    subplot(2,1,2)
    plot(sigmaP(2,:))
    
    figure
    subplot(2,1,1)
    plot(errpose(1,:))
    subplot(2,1,2)
    plot(errpose(2,:))
else disp('Done')
end