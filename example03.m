function example03
% this example runs an ekf on a setup with linear motion model 
% and nonlinear measurement model

% NOTE
% >> in the script means a big change is needed for you to adapt it to your
% setup. e.g. a new likelihood function or measurement model
% ^^ means there are parameters that can be tuned to see if it changes the
% performance

% generating 2D data from a mobile robot with unicycle kinematics
% set the parameters
dt=1/5;
tspan = 0:dt:20; % time for simulation

% ode with initial conditions 
i1=[2 2 0];
v=10; 
% this is just a simulation so we keep the unicycle model because it makes
% nice tracks.
[~,gt] = ode45(@unicycle_ode,tspan,i1,[],v);

% measurement model is for a radar
eta=1;
Z=radar1(gt(:,1:2)')+randn(2, size(gt,1))*eta;

% ready to filter
% initialize
X0(1:2,1)=gt(1,1:2); X0(3:4,1)=1; % the theta is just arbitrary
P0=eye(4)*.5; % ^^ initial covariance

try_extended_kalman_filter(dt, Z, gt', X0, P0);

function try_extended_kalman_filter(dt, Z, gt, X0, P0)

n=size(X0,1);
m=size(Z,1);
T=size(Z,2);


Xh=zeros(n,T); 
Xh_=Xh;
P=zeros(n,n,T); P_=P;


% kalman filter matrices
% >> constant velocity motion model and its ^^ process noise 
% Eq. 6.2.2-13 Bar-Shalom, RongLi, Kirubarajan, Estimation with
% Applications to Tracking and Navigation
Fk=  [ 1, 0, dt,  0
    0, 1, 0,  dt
    0, 0, 1,  0
    0, 0, 0,  1];

Ffun=@(x) Fk*x;

% ^^ this is where we plug in the Jacobian of a nonlinear motion model 
% evaluated at x  
Flinfun=@(x) Fk; 

Qk= eye(4)*1;

% >> measurement model
Hfun=@(x) radar1(x);

% ^^ this is where we plug in the Jacobian of a nonlinear measurement model
Hlinfun=@(x) radar1_lin(x); % because Hfun is linear

Rk=diag([2 2]*3);


k0=1;
kF=T;

% mean
Xh_(:,k0)=X0;

% covariance
P_(:,:,k0, 1)=P0;

% kalman filter
for k=k0:kF
    
    % update
    [Xh(:,k), P(:,:,k,1)]=ekf_update(Xh_(:,k), ...
        P_(:,:,k,1), Rk, Z(:,k), Hfun, Hlinfun);
    
    % predict
    [Xh_(:,k+1), P_(:,:,k+1,1)]= ekf_predict(Xh(:,k), ...
        P(:,:,k,1), Qk, Ffun, Flinfun);
end

show_the_results(gt, Z, Xh, P);

function show_the_results(gt, Z, Xh, P)

m=size(Z,1);

% show the results
figure(1); gcf; clf;
ylabels={'r_1', 'r_2'};
for ii=1:2
    subplot(2,1,ii); gca;
    plot(gt(ii,:), 'k', 'linewidth', 2);
    hold on;
%     plot(Z(ii,:), 'k*');
    plot(Xh(ii,:), 'r', 'linewidth', 2);
    set(gca, 'fontname', 'times', 'fontsize', 24);
    grid on;
    ylabel(ylabels{ii});
end
xlabel('time');


figure(2); gcf; clf;

plot(gt(1,:), gt(2,:), 'k', 'linewidth', 2);
hold on;
% plot(Z(1,:), Z(2,:), 'k*');
plot(Xh(1,:), Xh(2,:), 'r', 'linewidth', 2);
axis image;

set(gca, 'fontname', 'times', 'fontsize', 24);
legend('ground truth','estimate');


function Xdot = unicycle_ode(t,X,v)
Xdot(1,1) = v*cos(X(3));
Xdot(2,1) = v*sin(X(3));
% set omega as a function of time to create interesting trajectory
Xdot(3,1) = 1*sin(.5*t);

function z = radar1(X)

z(1,:)=sqrt(X(1,:).^2+X(2,:).^2);
z(2,:)=atan2(X(2,:), X(1,:));

function Hk = radar1_lin(X)

Hk=[X(1,1)/(sqrt(X(1,1)^2+X(2,1)^2)), X(2,1)/(sqrt(X(1,1)^2+X(2,1)^2)), 0 0
    1/(1+(X(2,1)/X(1,1))^2)*(-X(2,1)/X(1,1)^2), 1/(1+(X(2,1)/X(1,1))^2)*(1/X(1,1)) 0  0]; 


        
