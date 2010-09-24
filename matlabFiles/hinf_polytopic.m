%========================================
% Sistemas Multivariáveis - UFRGS
% Tassiano Neuhaus - tassianors@gmail.com
% Setembro 2010
%
% LMIs - Hinf - Polytopic
%========================================

% init
clear all;close all;

% limits
a=[-0.25 -2];
k=[-0.4649e-4 -0.7449e-4];
b=-0.012725;

% middle of uncertainties
A0=[[0 1];[-b*(a(1)+a(2))/2 b+(a(1)+a(2))/2]]
B0=[0;(k(1)+k(2))/2]

% System description - polytopic
A1=[[0 1];[-b*a(1) b+a(1)]]
A2=[[0 1];[-b*a(2) b+a(2)]]
B1=[0;k(1)]
B2=[0;k(2)]
Bw=B1;%for convinience
% I choose this value
C=[1 0];
 
% z=(G-HK)x
G=C;
H=0;
% get usefuul dimensions
[n,m]=size(B1);
%========================================
% Init LMI's descriptions
setlmis([]);

% var descriptions
W=lmivar(1,[n 1]);
R=lmivar(2,[m,n]);
gamma=lmivar(1,[1,1]);
%========================================
% equation 1 - A1 B1
lmi1a=newlmi;
% a11
lmiterm([-lmi1a 1 1 W],-A1,1,'s'); % Symmetric
lmiterm([-lmi1a 1 1 R],B1,1, 's'); % Symetric 
lmiterm([-lmi1a 1 1 0],-Bw*Bw');
% a12 = (G*W-H*R)'
lmiterm([-lmi1a 1 2 W],-1,G');      % W*G'
lmiterm([-lmi1a 1 2 -R],1,H');    % -R'*H'
%a22
lmiterm([-lmi1a 2 2 gamma],1,1);         % I
%========================================
% equation 2 - A1 B2
lmi1b=newlmi;
% a11
lmiterm([-lmi1b 1 1 W],-A1,1,'s'); % Symmetric
lmiterm([-lmi1b 1 1 R],B2,1, 's'); % Symetric 
lmiterm([-lmi1b 1 1 0],-Bw*Bw');
% a12 = (G*W-H*R)'
lmiterm([-lmi1b 1 2 W],-1,G');      % W*G'
lmiterm([-lmi1b 1 2 -R],1,H');    % -R'*H'
%a22
lmiterm([-lmi1b 2 2 gamma],1,1);         % I
%========================================
% equation 3 - A2 B1
lmi1c=newlmi;
% a11
lmiterm([-lmi1c 1 1 W],-A2,1,'s'); % Symmetric
lmiterm([-lmi1c 1 1 R],B1,1, 's'); % Symetric 
lmiterm([-lmi1c 1 1 0],-Bw*Bw');
% a12 = (G*W-H*R)'
lmiterm([-lmi1c 1 2 W],-1,G');      % W*G'
lmiterm([-lmi1c 1 2 -R],1,H');    % -R'*H'
%a22
lmiterm([-lmi1c 2 2 gamma],1,1);         % I
%========================================
% equation 4 - A2 B2
lmi1d=newlmi;
% a11
lmiterm([-lmi1d 1 1 W],-A2,1,'s'); % Symmetric
lmiterm([-lmi1d 1 1 R],B2,1, 's'); % Symetric 
lmiterm([-lmi1d 1 1 0],-Bw*Bw');
% a12 = (G*W-H*R)'
lmiterm([-lmi1d 1 2 W],-1,G');      % W*G'
lmiterm([-lmi1d 1 2 -R],1,H');    % -R'*H'
%a22
lmiterm([-lmi1d 2 2 gamma],1,1);         % I
%========================================

lmi2=newlmi;
lmiterm([-lmi2 1 1 W],1,1);      % W > 0
%========================================
% finish relationship description
nomesis = getlmis;

% solve
no = decnbr(nomesis);
co = zeros(no,1);
for j=1:no,
	[gammaj] = defcx(nomesis,j,gamma);
	co(j) = gammaj;
end

%Executa a minimiza�ao do criterio
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
A1_B2=ss(A1-B2*K,B2,C,0)
A2_B2=ss(A2-B2*K,B2,C,0)
A2_B1=ss(A2-B1*K,B1,C,0)
step(A0_B0,A1_B1, A1_B2, A2_B1,  A2_B2)
