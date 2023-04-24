clear
clc
global M1 M K b A1 B

%=== Variables ====================================================
M1=3.4;
%%%%%%%%%%%%
x0=[0.1 0 0 0 0 0];
ts=[0:2000];
%==================================================================

%=== Differential Eqn. ============================================
[t,x]=ode45('aero',ts,x0);
%==================================================================

%=== Result plot ==================================================
xi=x(:,1);
alpha=x(:,2);
%beta=x(:,3);
%dxi=x(:,4);
%dalpha=x(:,5);
%dbeta=x(:,6);

%plot(t,xi,'k:');



%******************************************************************

%--- SMC 2006.summer. ---------------------------------------------

%******************************************************************

[nn,mm]=size(B);

[Tr temp]=qr(B);                          % QR decomposition
Tr=Tr';
Tr=[Tr(mm+1:nn,:);Tr(1:mm,:)];

clear temp

%=== Regular form. =================================================

Areg=Tr*A1*Tr';
Breg=Tr*B;

%=== matrix sub-block. =============================================

A11=Areg(1:nn-mm,1:nn-mm);
A12=Areg(1:nn-mm,nn-mm+1:nn);
A21=Areg(nn-mm+1:nn,1:nn-mm);
A22=Areg(nn-mm+1:nn,nn-mm+1:nn);
B2=Breg(nn-mm+1:nn,1:mm);


%=== define the sliding Function using quadratic minimization. =====

% Positive definite Q
[Q,A12tilde]=qr(A12');

% Transform weighting matrix to regular form coordinate system
Qt=Tr*Q*Tr';   % 앞서 Tr을 변환해주었는데 이것 때문에 오류가 생긱수도 있음

% Compatibly partition with regular form description
Q11=Qt(1:nn-mm,1:nn-mm);
Q12=Qt(1:nn-mm,nn-mm+1:nn);
Q21=Qt(nn-mm+1:nn,1:nn-mm);
Q22=Qt(nn-mm+1:nn,nn-mm+1:nn);

% Form reduced order system description
Qhat=Q11-Q12*inv(Q22)*Q21;
Ahat=A11-A12*inv(Q22)*Q21;

% Solve the LQR problem
[K,P1,E]=lqr(Ahat,A12,Qhat,Q22);

% Obtain the switching function
M=inv(Q22)*(A12'*P1+Q21);
S2=eye(mm);             
S=S2*[M S2]*Tr;      
