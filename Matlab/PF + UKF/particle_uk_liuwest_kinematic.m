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

x = [0;0; 0; 0];
params = [1.5; 1.5];

%The system model for state UKF's
xsbar = [0; 0; 0; 0]; % initial belief
Ps = diag([1, 1, 0.2, 2]); %initial belief uncertainity
Rs = diag([0.5^2 0.5^2 0.2^2 1^2]); %process noise

Q = 100*diag([1^2 1^2 0.1^2 0.1^2]); %measurement noise
errpose = [];

%Parameter initialization for the particle filter
parameter_rangeH = [5; 5];
parameter_rangeL = [1; 1];
M = 1000; %number of particles

param_smooth = 0.009;
param_alpha = sqrt(1-param_smooth*param_smooth);
param_gamma = 0.05;

for i=1:length(parameter_rangeH)
    S.X(:,i) = (parameter_rangeH(i,1)-parameter_rangeL(i,1)).*rand(M,1) + parameter_rangeL(i,1);
end
S.W = (1/M)*ones(M,1); %initializing equal weights

X(:,1) = x;
Xsbar(:,1) = xsbar;
Xpbar(:,1) = S.X'*S.W; %weighted average

sim_noise = 0.1;
disp('Estimating parameters from kinematic model')
n = 1000; %number of iteration
for k = 1:n
    u = [2.5*cos(0.001*k),0.15*sin(0.001*k)];
    %u = [7.5e5,0.1];
    x = kinematic_step(x,params,u,dt,sim_noise);
    z = x + 0.05*randn(4,1); %observing all state in simulation mode
    diffusion = [];
    %% Parameter prediction using particle filter based on best state estimate available
    for i=1:M
        z_phat = kinematic_step(xsbar,S.X(i,:),u,dt,0);
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
        S.X(i,:) = m(i,:) + param_smooth*param_smooth*randn(1,2)*V;
    end
    xphat = S.W'*S.X
    %[~, idx] = sort(S.W, 'descend');
    %xphat = S.X(idx(1),:)
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
%     xphat = mean(S.X)
    %xphat = [xphat(1)*100,xphat(2)*1e4,xphat(3)*1e4]
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
        Ys(:,l)=kinematic_step(Xs(:,l),xphat,u,dt,0); %The map of state
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
    
    Xsbar(:,k+1) = xsbar;
    Xpbar(:,k+1) = xphat;
    X(:,k+1) = x;
    err = xphat - [1.5;1.5];
    errpose = [errpose err];
end
figure,
plot(Xpbar(1,:),'r')
figure,
plot(Xpbar(2,:),'b')
disp('Done')

figure
subplot(2,1,1)
plot(errpose(1,:))
subplot(2,1,2)
plot(errpose(2,:))
[Xpbar(1,n),Xpbar(2,n)]
[abs(Xpbar(1,n)-1.5)/1.5*100, abs(Xpbar(2,n)-1.5)/1.5*100]

