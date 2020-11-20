%% Dual UKF based state and parameter estimation for autonomous vehicle control
%
clc
clear all
close all


%% Estimation of Tire Co-efficients and Mass using Dynamic model
load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));
dt= 0.1;

pre_computed_params = [1.56,1.385];
known_params = [1.5,1.5];
%dynamic states = vx vy yaw_rate
%observing vx vy yaw_rate
x = [1;0;0];
%tire stiffness co-efficients
param = [1700; 15e4; 4e4];

%Parameter initialization for UKF
xsbar = [0.85;0;0]; % initial belief of state
Ps = diag([1, 1, 0.2]); %initial belief uncertainity
Rs = diag([1.2^2 1.2^2 0.45^2]); %process noise

Q = 10*diag([0.02^2 0.02^2 0.025^2]); %measurement noise
Hs = eye(3);
err_p = [];
err_s = [];

%Parameter initialization for the particle filter
parameter_rangeH =[20; 20; 8];
parameter_rangeL = [10; 10; 2];
M = 1000; %number of particles

param_smooth = 0.02;
param_alpha = sqrt(1-param_smooth*param_smooth);
param_gamma = 0.95;

for i=1:length(parameter_rangeH)
    S.X(:,i) = (parameter_rangeH(i,1)-parameter_rangeL(i,1)).*rand(M,1) + parameter_rangeL(i,1);
end
S.W = (1/M)*ones(M,1); %initializing equal weights

X(:,1) = x;
Xsbar(:,1) = xsbar;
Xpbar(:,1) = S.X'*S.W; %weighted average

sim_noise = 0.05;
disp('Estimating parameters from dynamic model')
n = size(u1,1);
for loop = 1:n
    loop;
    u = u1(loop,:);
    x = dynamic_step(x,u,param,known_params,dt,sim_noise);
    z = Hs*x + 0.05*randn(3,1);
    diffusion = 0;
    %% Parameter prediction using particle filter based on best state estimate available
    for i=1:M
        z_phat = dynamic_step(xsbar,u,[S.X(i,1)*100,S.X(i,2)*1e4,S.X(i,3)*1e4],pre_computed_params,dt,0);
        dz = z - z_phat;
        diffusion(i) = norm(z);
        S.W(i) = S.W(i)*exp(-0.5*dz'*inv(Q)*dz);
        if S.W(i) == 0
            S.W(i) = 1e-10;
        end
    end
    S.W = S.W/sum(S.W);
    xp_bar = S.W'*S.X;
    m = param_alpha*S.X + (1-param_alpha)*repmat(xp_bar,[M,1]);
    for i = 1:M
        V = S.W(i)*(S.X(i,:)-xp_bar)'*(S.X(i,:)-xp_bar);
        S.X(i,:) = m(i,:) + param_smooth*param_smooth*randn(1,3)*V;
    end
    xphat = S.W'*S.X;
    S.W = S.W*param_gamma;
    %%Systemaic sampling based on the weights
%     cdf = cumsum(S.W);
%     r_0 = rand/M;
%     for m = 1 : M
%         i = find(cdf >= r_0,1,'first');
%         S.X(m,:) = S.X(i,:);
%         r_0 = r_0 + 1/M;
%     end
%     S.W = (1/M)*ones(M,1);
%     xphat = mean(S.X);
    xphat = [xphat(1)*100,xphat(2)*1e4,xphat(3)*1e4];
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
    
    %% state prediction using unscented transform
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
    xsbar_prev = xsbar;
    xsbar=xshat+Ks*(z-zshat);                              %state update
    Ps=Pshat-Ks*Pzhat*Ks';                                %covariance update
    
    Xsbar(:,loop+1) = xsbar;
    Xpbar(:,loop+1) = xphat;
    X(:,loop+1) = x;
    
    s_err = xsbar - x;
    err_s = [err_s s_err];
    p_err = xphat - [1700;15e4;5e4];
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

Xpbar(:,loop+1)
mean(abs(err_s),2)
errPmean = mean(abs(err_p),2)
disp('Done')

