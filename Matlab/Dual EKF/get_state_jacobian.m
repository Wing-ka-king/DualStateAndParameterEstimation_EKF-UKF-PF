function J = get_state_jacobian(state,params, u,dt)

theta = state(3);
v = state(4);
u_2 = u(2);
lf = params(1);
lr = params(2);

J(1,:) = [ 1, 0, -dt*v*sin(theta + atan((lr*tan(u_2))/(lf + lr))),                       dt*cos(theta + atan((lr*tan(u_2))/(lf + lr)))];
J(2,:) = [ 0, 1,  dt*v*cos(theta + atan((lr*tan(u_2))/(lf + lr))),                       dt*sin(theta + atan((lr*tan(u_2))/(lf + lr)))];
J(3,:) = [ 0, 0,                                                1, (dt*tan(u_2))/((lf + lr)*((lr^2*tan(u_2)^2)/(lf + lr)^2 + 1)^(1/2))];
J(4,:) = [ 0, 0,                                                0,                                                                   1];

end