clear;close all;clc

% read the parameter matrix
eval('aero_para_deriv'); % parameters to design SMO

tspan=[0:300];
init=[Tr*[0.1;0.1;0.1;0;0;0]]';
options=[];

[t,z]=ode45('aero_smo_f',tspan,init,options,A11,A12,A21,A22,B2,S,M,Tr,ALnon,Knon);

% figure(1)
% subplot(221),plot(t,z(:,1))
% subplot(222),plot(t,z(:,2))
% subplot(223),plot(t,z(:,3))
% figure(2)
% subplot(221),plot(t,z(:,4))
% subplot(222),plot(t,z(:,5))
% subplot(223),plot(t,z(:,6))
% figure(3)
% subplot(221),plot(z(:,1),z(:,4))
% subplot(222),plot(z(:,2),z(:,5))
% subplot(223),plot(z(:,3),z(:,6))

x=(Tr'*z')'; % z=Tr*x

figure(4)
subplot(221),plot(t,x(:,1))
subplot(222),plot(t,x(:,2))
subplot(223),plot(t,x(:,3))
figure(5)
subplot(221),plot(t,x(:,4))
subplot(222),plot(t,x(:,5))
subplot(223),plot(t,x(:,6))
figure(6)
subplot(221),plot(x(:,1),x(:,4))
subplot(222),plot(x(:,2),x(:,5))
subplot(223),plot(x(:,3),x(:,6))

