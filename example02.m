function example02
% this example runs an ekf on a setup with nonlinear motion model 
% and linear measurement model

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
[~,gt] = ode45(@unicycle_ode,tspan,i1,[],v);
% ^^ add noise with intensity w0; gt is ground truth

% measurement model is for an overhead camera
eta=1; h=1; fl=1; px=320; py=240;
Z=overhead_cam(gt(:,1:2)',fl,  px, py,h)+randn(2, size(gt,1))*eta;

% ready to filter
% initialize
X0(1:2,1)=gt(1,1:2); X0(3,1)=0; % the theta is just arbitrary
P0=eye(3)*.5; % ^^ initial covariance

try_extended_kalman_filter(dt, Z, gt', X0, P0);

function try_extended_kalman_filter(dt, Z, gt, X0, P0)

n=size(X0,1);
m=size(Z,1);
T=size(Z,2);


Xh=zeros(n,T); 
Xh_=Xh;
P=zeros(n,n,T); P_=P;


% kalman filter matrices
% unicycle model discrete form
v=1; omega=1;
Ffun=@(x) unicycle1(x,v,omega, dt);

% ^^ this is where we plug in the Jacobian of a nonlinear motion model 
% evaluated at x  

Flinfun=@(x) unicycle_lin(x,v, dt); 

Qk= eye(3)*1;

% >> measurement model and ^^ noise
% this one is simply identity
h=1; fl=1; px=320; py=240;
Hfun=@(x) overhead_cam(x, fl, px,py, h);

% ^^ this is where we plug in the Jacobian of a nonlinear measurement model
Hk=[fl/h 0 0;
    0 fl/h 0];

Hlinfun=@(x) Hk; % because Hfun is linear

Rk=diag([2 2]*1);


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
ylabels={'x_1', 'x_2', '\theta'};
for ii=1:3
    subplot(3,1,ii); gca;
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

function X = unicycle1(X,v,omega, dt)
X(1,1) = X(1,1) + v*cos(X(3,1))*dt;
X(2,1) = X(2,1) + v*sin(X(3,1))*dt;
X(3,1) = X(3,1) + omega*dt;

function flin = unicycle_lin(X,v, dt)

flin=[  1 0 -v*sin(X(3,1)*dt);
        0 1 v*cos(X(3,1)*dt);
        0 0 1];

function z = overhead_cam(X, fl, px,py, h)

z(1,:)=X(1,:)/h*fl - px;
z(2,:)=X(2,:)/h*fl - py;


        
