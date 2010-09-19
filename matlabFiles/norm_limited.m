%Trabalho - Controle Robusto 
% LMIs - Norm-Bounded

% init
clear all;close all;

% limits
a=[-0.25 -2];
k=[-0.4649e-4 -0.7449e-4];
b=-0.012725;

% System description
A0=[[0 1];[-b*(a(1)+a(2))/2 b+(a(1)+a(2))/2]]
B0=[0;(k(1)+k(2))/2]
Ea=[[(a(1)-a(2))/2 0];[0 (a(1)-a(2))/2]]
Eb=[[0 0];[0 0]]
D=[[0 0];[1 1]]

% get usefuul dimensions
[n,m]=size(B0);
%[p,n]=size(C);


% Init LMI's descriptions
setlmis([]);

% var descriptions
W=lmivar(1,[n 1]);
R=lmivar(2,[m,n]);

% realtionship description
lmi1=newlmi;
% a11= -A0+B0*R-W*A0'+R'*B0'-D*D'
% variables are R and W
%lmiterm([lmi1 1 1 0],-A0);       % -A0
lmiterm([lmi1 1 1 R],B0,1);      % +B0*R
lmiterm([lmi1 1 1 W],-1,A0');    % -W*A0'
lmiterm([lmi1 1 1 -R],1,B0');    % +R'*B0'
lmiterm([lmi1 1 1 0],-D*D');     % -D*D' 
% a12 = -(W*Ea+R*Eb')
lmiterm([lmi1 2 1 W],-1,Ea);     % -W*Ea 
lmiterm([lmi1 2 1 R],-1,Eb');    % -R*Eb'
%a22
lmiterm([lmi1 2 2 0],-1);        % -I

lmi2=newlmi;
lmiterm([-lmi2 1 1 W],1,1);      % W > 0

% finish relationship description
nomesis = getlmis;

% solve
[aa,gopt]=feasp(nomesis);

% Get numerical results
W_o=dec2mat(nomesis,gopt,W)
R_o=dec2mat(nomesis,gopt,R)

% get K feedback matrix
K=R_o*inv(W_o)

% compare solutions
autoA=eig(A0)
autoABK=eig(A0+B0*K)