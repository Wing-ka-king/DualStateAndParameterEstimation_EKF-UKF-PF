function G = get_params_jacobian(state,params, u,dt)

theta = state(3);
v = state(4);
u_1 = u(1);
u_2 = u(2);
lf = params(1);
lr = params(2);

G(1,:) = [                                          (dt*lr*v*tan(u_2)*sin(theta + atan((lr*tan(u_2))/(lf + lr))))/((lf + lr)^2*((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)),                                                                             -(dt*v*sin(theta + atan((lr*tan(u_2))/(lf + lr)))*(tan(u_2)/(lf + lr) - (lr*tan(u_2))/(lf + lr)^2))/((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)];
G(2,:) = [                                         -(dt*lr*v*cos(theta + atan((lr*tan(u_2))/(lf + lr)))*tan(u_2))/((lf + lr)^2*((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)),                                                                              (dt*v*cos(theta + atan((lr*tan(u_2))/(lf + lr)))*(tan(u_2)/(lf + lr) - (lr*tan(u_2))/(lf + lr)^2))/((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)];
G(3,:) = [ (dt*lr^2*v*tan(u_2)^3)/((lf + lr)^4*((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)^(3/2)) - (dt*v*tan(u_2))/((lf + lr)^2*((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)^(1/2)), (dt*v*tan(u_2)*((2*lr^2*tan(u_2)^2)/(lf + lr)^3 - (2*lr*tan(u_2)^2)/(lf + lr)^2))/(2*(lf + lr)*((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)^(3/2)) - (dt*v*tan(u_2))/((lf + lr)^2*((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)^(1/2))];
G(4,:) = [                                                                                                                                                        0,                                                                                                                                                                                                                   0];

end