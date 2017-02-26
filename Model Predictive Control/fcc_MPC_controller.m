%------------ FCC: Model Predictive Controller and State Estimation -------
%
%              System  :   X(k+1) = A X(k)  +  B U(k) + [e1 ~ N(0,Q)]
%                          Y(k)   = C X(k)            + [e2 ~ N(0,R)] 
%
%                      
%
%     X-state         : (state-1 state-2 state-3 state-4)-[4 X 1]
%     U-control input : (F_a F_sc)                        [2 X 1]
%     Y- measurement  : Obtained from fccDynamics.m -     [3 X 1]
%     
%
%------------------------------------------------------------------------
% MATLAB R2014a

% Author : Koundinya 
%          AE13B010

fprintf('Loading plant Parameters\n');
load 'linssmodel.mat';

% Sampling Time
T_s=Ts;

% Run time-------- in terms of sampling time(T_s)
max_time_steps=400;

%---------System/Plant Model Parameters----
A=A;
B=B;
C=C_new;

% Number of states
  num_states=length(A);
% Number of control-inputs
  num_controlInputs=size(B,2);
% Number of measurements
  num_measurements=size(C,1);
  
  Y_measurement=zeros(num_measurements,max_time_steps);
%-----------MPC Controller Parameters--------

% Control Horizon
c=2;

% Prediction Horizon
p=8;

% Weights for cost function : no.of measurements X measurement horizon
weights=[0 0 120];

fprintf('Initializing controller and filter/n');
% Weighted matrix
weights_total=[];
for i=1:p
    weights_total=[weights_total weights/2^i];
end
phi=diag(weights_total);

% FCC plant Setpoint
Y_setPoint = [0.03 0.25 300];

% Initial State
X_initial=X0;

U_mpc=zeros(num_controlInputs,c);

%-----------Kalman Filter Parameters---------

% State Belief Mean 
bel_x=zeros(num_states,max_time_steps);

bel_x(:,1)=X0;

% State Belief Covariance
P=[1 0 0 0; 0 0 1 0 ;0 0 1 0;0 0 0 1];

% Measurement noise covariance
R = [0.9 0 0 ;0 0.1 0 ;0 0 1];

% Covariance of unmodelled random processes in system
Q = [0.1 0 0 0;0 0.5 0 0 ;0 0 0.6 0;0 0 0 1];

%------------------------------------------------
% vector to store result of optimization
vector_optimization=zeros(1,p);

% Array containing future states obtained from cost function optimization 
X_bar=zeros(length(X_initial),p);

k=1;

fprintf('Iterations started...\n');

while(k<=max_time_steps)
    fprintf('time step : %d\n',k);
    % calculate U @ t=K - given state estimate at t=k 
    % -> solve optimization prolem , find Yks
    f=[];
    for i=1:p
    f=[f Y_setPoint];
    end
    vector_optimization=quadprog(2*phi,transpose(-2*f*phi));
    
    % -> Find Xks from Yks
     for i=1:p
         Y_optimal=[];
         for j=1:num_measurements
             Y_optimal=[Y_optimal;vector_optimization(j)];
         end
     X_bar(:,i)=linsolve(C_new,Y_optimal);
     end
       
    % -> Find Uk from Xk+1 and Xk
     for i=1:c
         U_mpc(:,i)=linsolve(B,X_bar(:,i+1)-A*X_bar(:,i));
     end
    
    % pass input and simulate/sample measurement
    tspan=T_s*(k):(T_s*(k+1)-T_s*(k))/2 :T_s*(k+1);
    if k==1
        boundary_condition=C*X_initial;
    else
        boundary_condition=transpose(Y_measurement(:,k-1));
    end
    [t_ode,X]=ode45(@(tspan,x) fccDynamics(tspan,x,transpose(U_mpc(:,1))),tspan,boundary_condition);
     Y_measurement(:,k)= transpose(X(end,:));
     
     if(k~=max_time_steps)
    % obtain state estimate at t=K+1
    X_prediction=A*bel_x(:,k)+B*U_mpc(:,1);
    P=A*P*transpose(A)+Q;
    S=C*(P)*transpose(C)+R;
    K=P*transpose(C)/(S);
    bel_x(:,k+1)=X_prediction+K*(Y_measurement(:,k)-C*X_prediction);
    P=(eye(num_states)-K*C)*P;
     end
    
    k=k+1;
end
fprintf('Done \n');
t=1:1:max_time_steps;
%figure('name','State-3')
%plot(t,bel_x(3,:));
figure('name','MPC Controller Performance: Output 3')
plot(t,Y_measurement(3,:),'--r',t,Y_setPoint(3));
%figure('name','MPC Controller Performance: Output 2')
%plot(t,Y_measurement(2,:),'--r',t,Y_setPoint(2));
%figure('name','MPC Controller Performance: Output 1')
%plot(t,Y_measurement(1,:),'--r',t,Y_setPoint(1));
