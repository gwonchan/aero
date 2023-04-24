function [xdot,t]=aero_smo_f(t,xx1,options,A11,A12,A21,A22,B2,S,M,Tr,ALnon,Knon)

rho=10^-2; % un
Phi=-10^-3; % ul 
% extract the states
x=xx1(1:5);
x_res=xx1(6);
% special form of A tilde
A=[A11 A12;
   A21 A22];
% sliding function
S1=S(1,1:5);
S2=S(1,6);
P2=lyap(Phi,eye(1));
B=[zeros(5,1);B2];
% control input
ul=-inv(S*B)*S*A*xx1+inv(S*B)*Phi*S*xx1;
un=-rho*inv(S2*B2)*((P2*S*xx1)/norm(P2*S*xx1));
u=ul+un;
Trtran=Tr';
Knonz=Tr*Knon*Tr'*xx1*(Trtran(2,:)*xx1)^2;
ALnonz=Tr*ALnon*(Trtran(2,:)*xx1)^3;

xdot1=A11*x+A12*x_res+Knonz(1:5)+ALnonz(1:5);
xdot2=A21*x+A22*x_res+Knonz(6)+ALnonz(6)+B2*ul;

xdot = [xdot1; xdot2];
t
