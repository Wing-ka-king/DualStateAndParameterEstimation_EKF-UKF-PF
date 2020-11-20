function G = get_dynparam_jacobian(x, u, p, kp, dt)

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


G(1,:) = [                                                                                                              (dt*(tyreCoeff_y*sin(steerAngle)*(steerAngle - atan((vel_y + distCog_f*yawRate)/vel_x)) - tyreCoeff_x*cos(steerAngle)*(slipWheel_fl + slipWheel_fr)))/m^2,                                           (dt*cos(steerAngle)*(slipWheel_fl + slipWheel_fr))/m,                                                                                                -(dt*sin(steerAngle)*(steerAngle - atan((vel_y + distCog_f*yawRate)/vel_x)))/m];
G(2,:) = [                                                       -(dt*(tyreCoeff_x*sin(steerAngle)*(slipWheel_fl + slipWheel_fr) - tyreCoeff_y*atan((vel_y - distCog_r*yawRate)/vel_x) + tyreCoeff_y*cos(steerAngle)*(steerAngle - atan((vel_y + distCog_f*yawRate)/vel_x))))/m^2,                                           (dt*sin(steerAngle)*(slipWheel_fl + slipWheel_fr))/m,                                                    -(dt*(atan((vel_y - distCog_r*yawRate)/vel_x) - cos(steerAngle)*(steerAngle - atan((vel_y + distCog_f*yawRate)/vel_x))))/m];
G(3,:) = [ -(dt*(distCog_f*(tyreCoeff_x*sin(steerAngle)*(slipWheel_fl + slipWheel_fr) + tyreCoeff_y*cos(steerAngle)*(steerAngle - atan((vel_y + distCog_f*yawRate)/vel_x))) + distCog_r*tyreCoeff_y*atan((vel_y - distCog_r*yawRate)/vel_x)))/(m^2*(distCog_f/2 + distCog_r/2)^2), (distCog_f*dt*sin(steerAngle)*(slipWheel_fl + slipWheel_fr))/(m*(distCog_f/2 + distCog_r/2)^2), (dt*(distCog_r*atan((vel_y - distCog_r*yawRate)/vel_x) + distCog_f*cos(steerAngle)*(steerAngle - atan((vel_y + distCog_f*yawRate)/vel_x))))/(m*(distCog_f/2 + distCog_r/2)^2)];

end