clc
clear all
close all

syms x y theta v m lf lr u_1 u_2 dt

assume(lf,'Real')
assume(lr,'Real')
assume(x,'Real')
assume(y,'Real')
assume(theta,'Real')
assume(v,'Real')
assume(dt,'Real')
assume(u_1,'Real')
assume(u_2,'Real')

beta = atan((lr*tan(u_2))/(lf+lr));
eqn_1 = x + v*cos(theta+beta)*dt;
eqn_2 = y + v*sin(theta+beta)*dt;
eqn_3 = theta + (v*sin(beta)/lr)*dt;
eqn_4 = v + u_1*dt;

J(1,:) = jacobian(eqn_1,[x, y, theta, v]);
J(2,:) = jacobian(eqn_2,[x, y, theta, v]);
J(3,:) = jacobian(eqn_3,[x, y, theta, v]);
J(4,:) = jacobian(eqn_4,[x, y, theta, v]);


G(1,:) = jacobian(eqn_1,[lf, lr]);
G(2,:) = jacobian(eqn_2,[lf, lr]);
G(3,:) = jacobian(eqn_3,[lf, lr]);
G(4,:) = jacobian(eqn_4,[lf, lr]);
