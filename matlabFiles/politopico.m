%Trabalho - Controle Robusto 

% LMIs - Incertezas politópicas

% init
clear all;close all;

a=[-0.25 -2];
k=[-0.4649e-4 -0.7449e-4]
b=-0.012725;

A1=[[0 1];[-b*a(1) a(1)+b]];
A2=[[0 1];[-b*a(2) a(2)+b]];
B1=[0;k(1)];
B2=[0;k(2)];
C=[1 1]
s1 = ltisys(A1, B1, C, 0)
s2 = ltisys(A2, B2, C, 0);
polsys = psys([s1 s2])

% To test its quadratic stability, type
[tmin,P] = quadstab(polsys)

sys=ss(A1-B1*(P*B1)',B1,C,0)
step(sys)