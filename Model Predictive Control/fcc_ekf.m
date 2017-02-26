%------------ State Estimation of FCC using Extended Kalman Filter-------
%
%              System  :   X_dot = g(u,X) +  e1 (coavriance-Q)
%                          Y     = C h_a(X) +  e2 (covariance-R)
%
%                     h(X)=C*h_a(X)    g(u,X)-> given in fccDynamics.m   
%
%     X-state         : (C_rc  O_d  T_rg)
%     u-control input : (F_a F_sc)
%     Y- measurement 
%
%------------------------------------------------------------------------
% MATLAB R2014a

% Author : Koundinya 
%          AE13B010


%----Initializing system model variables----


C_1 = [0 0 1];           % one measured output 
C_2 = [1 0 0; 0 0 1];    % two measured outputs

syms x1 x2 x3
h_sym_case1 = C_1*[x1;x2;x3];
h_sym_case2 = C_2*[x1;x2;x3];

%--- Initializing Kalman Filter variables---
bel_x=zeros(3,4094);
P_t=zeros(3,3);
P_T=zeros(3,3);

%R=zeros(2,4094);
%Q=zeros(3,4094);
%X=zeros(3,4);
%p=zeros(3,4);
K=zeros(3,2);
% Initial state belief mean
bel_x(1,1) = 0.035;
bel_x(2,1) = 0.001;
bel_x(3,1) = 900;

% Initial state belief covariance
 P_t=[1 0 0; 0 1 0 ;0 0 1];

% Measurement noise covariance
R = [0.9 0 ;0 0.1 ];

% Covariance of unmodelled random processes in system
Q = [0.1 0 0;0 0.5 0 ;0 0 0.6];

% maximum number of time-steps/measurement samples : 
 dt_max = 4094;
 
 t_ode=0;

 % time-step counter
 t=1;
 
 Time = [Time ;-1];
 
 while(t<4095)
     fprintf('--Step :%d--\n',t);
 % EKF Prediction
 tspan=Time(t):(Time(t+1)-Time(t))/2 :Time(t+1);
 
 [t_ode,X]=ode45(@(tspan,x) fccDynamics(tspan,x,InputU(t,:)),tspan,transpose(bel_x(:,t)));
 %[t_ode,p]=ode45(,tspan,P(:,t));
 bel_x(:,t)= transpose(X(end,:));
 fprintf('ODE solved \n');
 
 syms x1 x2 x3 u1 u2
 F_sym = jacobian (fccDynamics(t,[x1 x2 x3],[u1 u2]),[x1,x2,x3]);
 H_sym = jacobian (h_sym_case2,[x1 x2 x3]);
 
 F     =  subs(F_sym,[x1 x2 x3 u1 u2],[bel_x(1,t) bel_x(2,t) bel_x(2,t) InputU(t,1) InputU(t,2)]);
 H     =  subs(H_sym,[x1 x2 x3],[bel_x(1,t) bel_x(2,t) bel_x(3,t)]);
 fprintf('Jacobian found\n');
 P_T = P_t + (F*P_t + P_t*transpose(F) + Q)*(Time(t+1)-Time(t));

 % EKF Update

 K     = (P_T*transpose(H))*inv(H*P_T*transpose(H)+R);
 fprintf('Kalman Gain Calculated\n');
 bel_x(:,t+1)= bel_x(:,t) + K*(Y_measured_case2(t,:).'- C_2*bel_x(:,t));
 P     = (eye(3,3)-K*H)*P_T;
 P_t   =   P_T;
 fprintf('State Update Done\n');
 t=t+1;
 
 end
Time(end)=[];
plot(Time,transpose(bel_x(1,:)),Time,Y_measured_case2(:,1))


