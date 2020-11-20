clc
clear all
close all

syms dt
syms m distCog_f distCog_r tyreCoeff_x tyreCoeff_y
assume(m,'Real')
assume(distCog_f,'Real')
assume(distCog_r,'Real')
assume(tyreCoeff_x,'Real')
assume(tyreCoeff_y,'Real')

Jz = m*((distCog_f+distCog_r)/2)^2;

syms vel_x vel_y yawRate
assume(vel_x,'Real')
assume(vel_y,'Real')
assume(yawRate,'Real')
syms slipWheel_fl slipWheel_fr slipWheel_rl slipWheel_rr steerAngle

slipAngle_f = steerAngle - atan((vel_y + distCog_f*yawRate)/vel_x);
slipAngle_r = -atan((vel_y - distCog_r*yawRate)/vel_x);

Fy_f = tyreCoeff_y*slipAngle_f;
Fy_r = tyreCoeff_y*slipAngle_r;
Fx_f = tyreCoeff_x*(slipWheel_fl+slipWheel_fr);
Fx_r = 0;

torqueLong = (Fy_f*cos(steerAngle) + Fx_f*sin(steerAngle))*distCog_f - Fy_r*distCog_r;
acc_x = (Fx_f*cos(steerAngle) - Fy_f*sin(steerAngle) + Fx_r)/m;
acc_y = (Fy_f*cos(steerAngle) + Fx_f*sin(steerAngle) + Fy_r)/m;

dvX_dt = (acc_x + vel_y*yawRate);
dvY_dt = (acc_y - vel_x*yawRate);
dyawRate_dt = (torqueLong/Jz);


eqn1 = vel_x + dvX_dt*dt;
eqn2 = vel_y + dvY_dt*dt;
eqn3 = yawRate + dyawRate_dt*dt;


J(1,:) = jacobian(eqn1,[vel_x, vel_y, yawRate]);
J(2,:) = jacobian(eqn2,[vel_x, vel_y, yawRate]);
J(3,:) = jacobian(eqn3,[vel_x, vel_y, yawRate]);

G(1,:) = jacobian(eqn1,[m, tyreCoeff_x, tyreCoeff_y]);
G(2,:) = jacobian(eqn2,[m, tyreCoeff_x, tyreCoeff_y]);
G(3,:) = jacobian(eqn3,[m, tyreCoeff_x, tyreCoeff_y]);