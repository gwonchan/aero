function aero=f(t,x)
load A11.mat
load A12.mat
load MM.mat
aero=(A11-A12*MM)*x;

