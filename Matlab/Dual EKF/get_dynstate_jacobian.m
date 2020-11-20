function J = get_dynstate_jacobian(x, u, p, kp, dt)

m = p(1);                 %       m  = p(1);   /* Vehicle mass.                    */
distCog_f = kp(1);         %       a  = p(2);   /* Distance from front axle to COG. */
distCog_r = kp(2);         %       b  = p(3);   /* Distance from rear axle to COG.  */
tyreCoeff_x = p(2);       %       Cx = p(4);   /* Longitudinal tire stiffness.     */
tyreCoeff_y = p(3);       %       Cy = p(5);   /* Lateral tire stiffness.          */

vel_x = x(1);             %       /* x(1): Longitudinal vehicle velocity. */
vel_y = x(2);             %       /* x(2): Lateral vehicle velocity. */
yawRate = x(3);           %       /* x(3): Yaw rate. */

slipWheel_fl = u(1);      %       u1(t) = s_FL(t)     Slip of Front Left tire (ratio).
slipWheel_fr = u(2);      %       u2(t) = s_FR(t)     Slip of Front Right tire (ratio).
slipWheel_rl = u(3);      %       u3(t) = s_RL(t)     Slip of Rear Left tire (ratio).
slipWheel_rr = u(4);      %       u4(t) = s_RR(t)     Slip of Rear Right tire (ratio).
steerAngle = u(5);        %       u5(t) = delta(t)    Steering angle (rad)

J(1,:) = [                                                                                                                                                        1 - (dt*tyreCoeff_y*sin(steerAngle)*(vel_y + distCog_f*yawRate))/(m*vel_x^2*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1)),                                                                                                                 dt*(yawRate + (tyreCoeff_y*sin(steerAngle))/(m*vel_x*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1))),                                                                                                                 dt*(vel_y + (distCog_f*tyreCoeff_y*sin(steerAngle))/(m*vel_x*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1)))];
J(2,:) = [                                           -dt*(yawRate - ((tyreCoeff_y*(vel_y - distCog_r*yawRate))/(vel_x^2*((vel_y - distCog_r*yawRate)^2/vel_x^2 + 1)) + (tyreCoeff_y*cos(steerAngle)*(vel_y + distCog_f*yawRate))/(vel_x^2*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1)))/m),                                                   1 - (dt*(tyreCoeff_y/(vel_x*((vel_y - distCog_r*yawRate)^2/vel_x^2 + 1)) + (tyreCoeff_y*cos(steerAngle))/(vel_x*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1))))/m,                                -dt*(vel_x - ((distCog_r*tyreCoeff_y)/(vel_x*((vel_y - distCog_r*yawRate)^2/vel_x^2 + 1)) - (distCog_f*tyreCoeff_y*cos(steerAngle))/(vel_x*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1)))/m)];
J(3,:) = [ -(dt*((distCog_r*tyreCoeff_y*(vel_y - distCog_r*yawRate))/(vel_x^2*((vel_y - distCog_r*yawRate)^2/vel_x^2 + 1)) - (distCog_f*tyreCoeff_y*cos(steerAngle)*(vel_y + distCog_f*yawRate))/(vel_x^2*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1))))/(m*(distCog_f/2 + distCog_r/2)^2), (dt*((distCog_r*tyreCoeff_y)/(vel_x*((vel_y - distCog_r*yawRate)^2/vel_x^2 + 1)) - (distCog_f*tyreCoeff_y*cos(steerAngle))/(vel_x*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1))))/(m*(distCog_f/2 + distCog_r/2)^2), 1 - (dt*((distCog_r^2*tyreCoeff_y)/(vel_x*((vel_y - distCog_r*yawRate)^2/vel_x^2 + 1)) + (distCog_f^2*tyreCoeff_y*cos(steerAngle))/(vel_x*((vel_y + distCog_f*yawRate)^2/vel_x^2 + 1))))/(m*(distCog_f/2 + distCog_r/2)^2)];
 
end