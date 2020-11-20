%   /* State equations. */
function x_next = dynamic_step(x, u, p, kp, dt, noise_level)
%       /* Retrieve model parameters. */
%       /* Retrieve model parameters. */
m = p(1);                 %       m  = p(1);   /* Vehicle mass.                    */
distCog_f = kp(1);         %       a  = p(2);   /* Distance from front axle to COG. */
distCog_r = kp(2);         %       b  = p(3);   /* Distance from rear axle to COG.  */
tyreCoeff_x = p(2);       %       Cx = p(4);   /* Longitudinal tire stiffness.     */
tyreCoeff_y = p(3);       %       Cy = p(5);   /* Lateral tire stiffness.          */
Jz = m*((distCog_f+distCog_r)/2)^2;

vel_x = x(1);             %       /* x(1): Longitudinal vehicle velocity. */
vel_y = x(2);             %       /* x(2): Lateral vehicle velocity. */
yawRate = x(3);           %       /* x(3): Yaw rate. */

slipWheel_fl = u(1);      %       u1(t) = s_FL(t)     Slip of Front Left tire (ratio).
slipWheel_fr = u(2);      %       u2(t) = s_FR(t)     Slip of Front Right tire (ratio).
slipWheel_rl = u(3);      %       u3(t) = s_RL(t)     Slip of Rear Left tire (ratio).
slipWheel_rr = u(4);      %       u4(t) = s_RR(t)     Slip of Rear Right tire (ratio).
steerAngle = u(5);        %       u5(t) = delta(t)    Steering angle (rad)

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

x_next = x + [dvX_dt+noise_level*randn(1,1); dvY_dt+noise_level*randn(1,1); dyawRate_dt+noise_level*randn(1,1)*pi/180]*dt;

end