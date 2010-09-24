%========================================
% Sistemas Multivariáveis - UFRGS
% Tassiano Neuhaus - tassianors@gmail.com
% Setembro 2010
%
% LMIs - Hinf - Norm-Bounded
%========================================

% init
clear all; close all;

% limits
a=[-0.25 -2];
k=[-0.4649e-4 -0.7449e-4];
b=-0.012725;

% System description
A0=[[0 1];[-b*(a(1)+a(2))/2 b+(a(1)+a(2))/2]]
B0=[0;(k(1)+k(2))/2]
Bw=B0;%for convinience
% I choose this value
C=[1 0];
Ea=[[(a(1)-a(2))/2 0];[0 (a(1)-a(2))/2]]
Eb=[0; (k(1)-k(2))/2]
D=[[0 0];[1 1]]

% System description - polytopic
A1=[[0 1];[-b*a(1) b+a(1)]]
A2=[[0 1];[-b*a(2) b+a(2)]]
B1=[0;k(1)]
B2=[0;k(2)]

% z=(G-HK)x
G=C;
H=0;
[n,m]=size(B0);

% Init LMI's descriptions
setlmis([]);

% var descriptions
W=lmivar(1,[n 1]);
R=lmivar(2,[m,n]);
gamma=lmivar(1,[1,1]);

% realtionship description
lmi1=newlmi;

% variables are R, W and gamma
lmiterm([-lmi1 1 1 W],-A0,1,'s'); % symmetric
lmiterm([-lmi1 1 1 R],B0,1,'s'); % symmetric
lmiterm([-lmi1 1 1 0],-Bw*Bw');   % Choose Bw = B0 for convenience (not specified) 
%lmiterm([-lmi1 1 1 0],-D*D');     % -D*D'
% a12 = (Ea*W-Eb*R)'
lmiterm([-lmi1 1 2 W],1,Ea');     % W*Ea'
lmiterm([-lmi1 1 2 -R],-1,Eb');   % -R'*Eb'
% a13 = (G*W-H*R)'
lmiterm([-lmi1 1 3 W],1,-G');      % W*G'
lmiterm([-lmi1 1 3 -R],1,H');    % -R'*H'
%a22
lmiterm([-lmi1 2 2 0],1);         % I
% A33
lmiterm([-lmi1 3 3 gamma],1,1);   % I

lmi2=newlmi;
lmiterm([-lmi2 1 1 W],1,1);      % W > 0

% finish relationship description
nomesis = getlmis;

% solve
no = decnbr(nomesis);
co = zeros(no,1);
for j=1:no,
	[gammaj] = defcx(nomesis,j,gamma);
	co(j) = gammaj;
end

% Execute criteria
[copt,gopt] = mincx(nomesis,co);

% Get numerical results
W_o=dec2mat(nomesis,gopt,W)
R_o=dec2mat(nomesis,gopt,R)
gamma_o=dec2mat(nomesis,gopt,gamma)

% get K feedback matrix
K=R_o*inv(W_o)

% compare solutions
autoA=eig(A0)
autoABK=eig(A0-B0*K)

A0_B0=ss(A0-B0*K,B0,C,0)
A1_B1=ss(A1-B1*K,B1,C,0)
A2_B2=ss(A2-B2*K,B2,C,0)
step(A0_B0,A1_B1, A2_B2)

%sys=ss(A0-B0*K,B0,C,0);
%step(sys);

