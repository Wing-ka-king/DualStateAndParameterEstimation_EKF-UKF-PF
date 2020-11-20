% Dual EKF based state and parameter estimation for autonomous vehicle control 

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

%The system model for both EKF's
xsbar = [0; 0; 0; 0]; % initial belief
xpbar = [4; 4]; % initial belief

Ps = diag([1, 1, 0.2, 2]); %initial belief uncertainity
Pp = diag([5, 5]);
Rs = diag([0.5^2 0.5^2 0.2^2 1^2]); %process noise
Rp = diag([0.5^2 0.5^2]);
Q = 10*diag([1^2 1^2 0.1^2 0.1^2]); %measurement noise

n = 1000; %number of iteration
Hs = eye(4);

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
    u = [2.5*cos(0.001*k),0.15*sin(0.001*k)];
    %u = [7.5e5,0.1];
    x = kinematic_step(x,params,u,dt,sim_noise);
    z = x + 0.1*randn(4,1); %observing all state in simulation mode
    
    %parameter prediction
    xphat = xpbar;
    Pp = Pp + Rp;
    
    %state prediction
    xshat = kinematic_step(xsbar,xphat,u,dt,0);
    
    G = get_state_jacobian(xsbar,xphat,u,dt);
    Ps = G*Ps*G' + Rs;
    
    %calculating Kalman gain for state 
    Ks = Ps*Hs'*inv(Hs*Ps*Hs' + Q);
    
    %update step for state
    xsbar = xshat + Ks*(z - Hs*xshat);
    Ps = Ps - Ks*Hs*Ps;
    
    %calculating Kalman gain for parameters
    Hp = get_params_jacobian(xsbar,xphat,u,dt);
    Kp = Pp*Hp'*inv(Hp*Pp*Hp' + Q);
    
    %update step for parameters
    xpbar = xphat + Kp*(z - Hs*xshat);
    Pp = Pp - Kp*Hp*Pp;
    
    Xsbar(:,k+1) = xsbar;
    Xpbar(:,k+1) = xpbar;
    X(:,k+1) = x;
    sigmaP = [sigmaP Pp(:)];
    sigmaS = [sigmaS Ps(:)];
    err = xpbar - [1.5;1.5];
    errpose = [errpose err];
end
% figure, 
% plot(Xsbar(1,:),'b')
% hold on;
% plot(X(1,:),'r')
% 
% figure,
% plot(Xsbar(2,:),'b')
% hold on;
% plot(X(2,:),'r')
% 
% 
% figure,
% plot(Xsbar(3,:),'b')
% hold on;
% plot(X(3,:),'r')

% 
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

[Xpbar(1,n), Xpbar(2,n)],[abs(Xpbar(1,n)-1.5)/1.5*100, abs(Xpbar(2,n)-1.5)/1.5*100]
disp('Done')