%clearing the workspace, command window and previous plots
clear all
close all
clc

%Assigning symbolic variables 
syms x1 x2; % defining the state variables
u = 1; % the input variable chosen arbitrarily
u0 = 2; k1=0.06; k2=0.3; k =2; a=0.7; %given Variables
g= ((u0*x2*x1)/(k1+x2+k2*(x2^2)))-x1*u %x1 dot
h= (-(u0*x2*x1)/(a*(k1+x2+k2*(x2^2))))+(k-x2)*u; %x2 dot
[solve_x1,solve_x2]=solve([g==0,h==0],[x1,x2],"real",true);


x1_roots=solve_x1; 
x2_roots=solve_x2;
x1_value=round(x1_roots(end),4) % rounding the values of x1 and x2 to the nearest 4th decimal number
x2_value=round(x2_roots(end),4) 
disp(x1_value,x2_value);

g = ((u0*x2*x1)/(k1+x2+k2*(x2^2)))-(x1*u);
h = ((-u0*x2*x1)/(a*(k1+x2+k2*(x2^2)))) + u*(k-x2);
der1=diff(g,x1);
der2=diff(g,x2);
der3=diff(h,x1);
der4=diff(h,x2);
syms x1 x2
x1=-0.89054;
x2=3.272;
u=1;
d1=round(vpa((2*x2)/((3*x2^2)/10 + x2 + 3/50) - u),4);
d2=round(vpa((2*x1)/((3*x2^2)/10 + x2 + 3/50) - (2*x1*x2*((3*x2)/5 + 1))/((3*x2^2)/10 + x2 + 3/50)^2 ),4);
d3=round(vpa(-(2*x2)/(a*((3*x2^2)/10 + x2 + 3/50))),4);
d4=round(vpa((2*x1*x2*((3*x2)/5 + 1))/(a*((3*x2^2)/10 + x2 + 3/50)^2) - (2*x1)/(a*((3*x2^2)/10 + x2 + 3/50)) - u),4);

%Forming a jacobian matrix from the partial derivaties with assigned operation point 
%Forming matrices according to the equations
%xdot=Ax+Bu
%y=Cx+Du
A=[d1 d2;d3 d4]
B=[-x1;k-x2]
C=[1 0] %Matrix C is a given matrix
D=[0] %Matrix D is a given matrix
A=[0 0.13109;-1.42861 -1.18727]
B=[  0.8905 ;  -1.2720]

%Assigning time
t=0:0.01:1;

%Forming a linear state space model
sys=ss(A,B,C,D);

%Finding out step response
figure
step(sys,t)

%Checking stability

r=rlocus(sys); % to check the poles of the system on the graph
pole(sys) % to see the exact numerical values of the poles

Stability=eig(A)
if Stability<0 
    disp('System is stable')
else
    disp('System is unstable')
end

%Checking observability
Observability=obsv(A,C); 
if (length(A) - rank(Observability)==0)% checking whether the rank of the observability is equal to the length of A
    
    disp('System is observable')
else
    disp('System is unobservable')
end

%Checking controllability
Controllability=ctrb(A,B);
if (length(A) - rank(Controllability)==0) % checking whether the rank of the controllabelity is equal to the length of A
    
      disp('System is controllable')
else
    disp('System is uncontrollable')
end
unco = length(A) - rank(Controllability);
display(unco); % double checking to see if there is any uncontrollable state

%Designing a state space controller using linear quadratic approach
Q=C'*C;
R=0.001;
nx=2;
ny=1;
[K, S, E]=lqr(A, B, Q, R)
Ku=K(1);
Kv=K(2);
Acl=A-B*K;
sysss=ss(Acl,B,C,D);
step(sysss,t);
Nbar=-inv(C*inv(Acl)*B);
BB=Nbar*B;
sys_cl=ss(Acl,BB,C,D);
figure();
step(sys_cl,t)
S1=stepinfo(sys_cl)

%Designing a state observer
poles_con=E;
poles_obs=[-300,-300];

%Pole placement for observer
L=acker(A',C',poles_obs)';

%Finding transfer function of full order observer-controller
A_bar=A-L*C-B*K;
B_bar=L;
C_bar=K;
D_bar=zeros(size(C_bar,1),size(B_bar,2));

[num,den]=ss2tf(A_bar,B_bar,C_bar,D_bar);
sys1=tf(num,den);

%Finding transfer function of original system
[num,den]=ss2tf(A,B,C,D);
sys2=tf(num,den);

%Transfer function of system with observer-controller
sys3=sys1*sys2;
sys4=feedback(sys3,1);
figure();
step(sys4,t);
figure();
h=stepplot(sys4,sys_cl);
S2=stepinfo(sys4)


