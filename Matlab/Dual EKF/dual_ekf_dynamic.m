%% Estimation of Mass and Tire Co-efficients and Mass using Dynamic model
clc
close all
clear all
%loading control input data
load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

dt= 0.1;
pre_computed_params = [1.56,1.35];
orig_params = [1.5,1.5];
%dynamic states = vx vy yaw_rate 
%observing vx yaw_rate
x = [1;0;0.2];
%tire stiffness co-efficients
param = [1700; 15e4; 5e4];

% The system model for both UKF's
xsbar = [0.8;0;0]; % initial belief of state
xpbar = [1900; 13e4; 4e4]; % initial belief of parameter

Ps = diag([1^2, 1^2, 0.5^2]); %initial belief uncertainity
Pp = diag([300 3e4 1e4]);
Rs = diag([1.2^2 1.2^2 0.45^2]); %process noise
Rp = diag([10^2 100^2 75^2]);

Q = 10*diag([0.02^2 0.02^2 0.025^2]); %measurement noise
Hs = eye(3);
gamma = 1;
n = size(u1,1);
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;
X(:,1) = x;
Ppp=[Pp(1,1);Pp(2,2)];

disp('Estimating parameters from dynamic model')
time = [];
current_time = 0;

for loop = 1:n
    loop;
    time = [time;current_time];
    current_time = current_time+dt;
    u = u1(loop,:);
    x = dynamic_step(x,u,param,orig_params,dt,0.001);
    z = Hs*x + 0.01*randn(3,1); 
    
%     parameter prediction
    xphat = xpbar;
    Pp = Pp + Rp;
    Pp = (Pp + Pp')/2;
%     state prediction
    xshat = dynamic_step(xsbar,u,xphat,pre_computed_params,dt,0);
   
    xshat-x;
    
    G = get_dynstate_jacobian(xsbar,u,xphat,pre_computed_params,dt);
    Ps = G*Ps*G' + Rs;
    Ps = (Ps + Ps')/2;
%     calculating Kalman gain for state 
    Ks = Ps*Hs'*inv(Hs*Ps*Hs' + Q);
    
%     update step for state
    xsbar = xshat + Ks*(z - Hs*xshat);
    Ps = Ps - Ks*Hs*Ps;
    
%     calculating Kalman gain for parameters
    Hp = get_dynparam_jacobian(xsbar,u,xphat,pre_computed_params,dt);
    
    Kp = Pp*Hp'*inv(Hp*Pp*Hp' + Q);
    
%     update step for parameters
    xpbar = xphat + Kp*(z - Hs*xshat);
    Pp = Pp - Kp*Hp*Pp;
    
    Xsbar(:,loop+1) = xsbar;
    Xpbar(:,loop+1) = xpbar;
    X(:,loop+1) = x;
 end
figure,
subplot(3,1,1)
plot(Xsbar(1,:),'r','LineWidth',1)
hold on
plot(X(1,:),'b')
subplot(3,1,2)
plot(Xsbar(2,:),'r','LineWidth',1)
hold on
plot(X(2,:),'b')
subplot(3,1,3)
plot(Xsbar(3,:),'r','LineWidth',1)
hold on
plot(X(3,:),'b')

figure
subplot(3,1,1)
plot(Xpbar(1,:))
subplot(3,1,2)
plot(Xpbar(2,:))
subplot(3,1,3)
plot(Xpbar(3,:))

Xpbar(:,loop+1)

disp('Done')
% disp('Done')
% 
% plot(time(1:end-28),Xsbar(3,30:end),'b')
% hold on;
% plot(time(1:end-28),X(3,30:end),'r')