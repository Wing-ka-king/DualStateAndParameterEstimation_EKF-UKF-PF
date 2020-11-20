function pred_sim = kinematic_step(state,params,u,dt,noise_level)
%params = [lf,lr]
%states = [x,y,phi,v];
%u = [acc,steering_angle]
%equations based on the paper Kong et al, Kinematic and Dynamic Vehicle
%Models for Autonomous Driving control design

pred_sim = zeros(size(state));
beta = atan((params(2)*tan(u(2)))/(params(1)+params(2)));
pred_sim(1,1) = state(1) + state(4)*cos(state(3)+beta)*dt + noise_level*randn(1,1);
pred_sim(2,1) = state(2) + state(4)*sin(state(3)+beta)*dt + noise_level*randn(1,1);
pred_sim(3,1) = state(3) + state(4)*sin(beta)/params(2)*dt + noise_level*randn(1,1)*pi/180;
%pred_sim(3,1) = mod(pred_sim(3,1)+pi,2*pi)-pi; %this is not working yet
pred_sim(4,1) = state(4) + u(1)*dt + noise_level*randn(1,1)*2;
end