%clear all
%clc
%% Input IMU Data
dir = 'C:\Users\dschul\Dropbox (University of Michigan)\ENGIN-PREACT\2 - Simulation\Monica Jones Dataset\DS_Copy\UMTRI Mcity Data\';
subjects_mcity = ['MS019_FN'];
subject_folder = subjects_mcity(1,:);
    
subject = subject_folder(1:5);
condition = subject_folder(7:8);

dataDir = [dir,subject_folder];
%[h,v] = Data_handle.readtsv(dataDir);
filename_h = strcat(subject_folder,'_h','.mat');
filename_v = strcat(subject_folder,'_v','.mat');
    
load(filename_h,'h')
load(filename_v,'v')

[ax,ay,az,wx,wy,wz,timeStop] = Data_handle.setVehAcc(v);

%% Correct for shift in IMU coordinate frame
k = [mean(ax); mean(ay); mean(az)];
g = [0;0;-9.8];

k_norm = k/norm(k);
g_norm = g/norm(g);
angle = acos(dot( k_norm, g_norm ));
axis = cross( k_norm, g_norm );
axis = axis/norm(axis);

c = cos(angle);
s = sin(angle);
t = 1 - c;
x = axis(1);
y = axis(2);
z = axis(3);
R = [[c+x*x*t , x*y*t-z*s, x*z*t+y*s]; [y*x*t+z*s, c+y*y*t, y*z*t-x*s]; [z*x*t-y*s, z*y*t+x*s, c+z*z*t]];  

adata = [ax.Data(:)'; ay.Data(:)'; az.Data(:)'];
adata = R*adata;

ax.Data = adata(1,:);
ay.Data = adata(2,:);
az.Data = -adata(3,:);

wx = wx - mean(wx);
wy = wy - mean(wy);
wz = wz - mean(wz);
    
%% Calculate desired trajectory - alpha and gamma to motor angles
%These values enforce experienced accel coincident with passenger y axis (sagital)

alpha = atan(ay.Data./9.8);
gamma = atan(ax.Data./9.8);
theta1_d = zeros(1,length(alpha)); %desired motor angle 1 (r1)
theta2_d = zeros(1,length(alpha)); % desired motor angle 2 (r2)

%Find motor angles based on inverse kinematic equations

G = [0 0 0];                %Origin
M = [0.23 0.29 0.07];       %M Measurements [meter]
S = [0.23 0.37 0];          %S Measurements [meter]
l_gl = 0;
l_mp = 0.06;
l_ps = 0.07;
m = 1;

for i = 1:length(alpha)

mu1 = S(1)*cos(alpha(i)) + S(3)*sin(alpha(i)) - M(1);
mu2 = -S(1)*cos(alpha(i)) + S(3)*sin(alpha(i)) + M(1);

sigma1 = S(1)*sin(alpha(i))*sin(gamma(i)) + S(2)*cos(gamma(i)) - S(3)*cos(alpha(i))*sin(gamma(i)) - M(2);
sigma2 = -S(1)*sin(alpha(i))*cos(gamma(i)) + S(2)*sin(gamma(i)) + S(3)*cos(alpha(i))*cos(gamma(i)) + M(3);
sigma3 = -S(1)*sin(alpha(i))*sin(gamma(i)) + S(2)*cos(gamma(i)) - S(3)*cos(alpha(i))*sin(gamma(i)) - M(2);
sigma4 = S(1)*sin(alpha(i))*cos(gamma(i)) + S(2)*sin(gamma(i)) + S(3)*cos(alpha(i))*cos(gamma(i)) + M(3);

a1 = sigma1;
a2 = sigma3;

b1 = sigma2;
b2 = sigma4;

c1 = ((l_ps^2) - (mu1^2) - (sigma1^2) - (sigma2^2) - (l_mp^2))/(-2*l_mp);
c2 = ((l_ps^2) - (mu2^2) - (sigma3^2) - (sigma4^2) - (l_mp^2))/(-2*l_mp);

phi1 = atan2(a1,b1);
phi2 = atan2(a2,b2);

r1 = sqrt((a1^2) + (b1^2));
r2 = sqrt((a2^2) + (b2^2));

theta1_d(i) = real(asin(c1/r1) - phi1);
theta2_d(i) = real(asin(c2/r2) - phi2);
end

%% Define continuous system dynamics (command voltage to motor angle)
J = 8.5e-6;
Bm = 4.2e-6;
Kt = 0.022945;
Ka = 2;

A = [0 1 0 0; 0 -Bm/J 0 0; 0 0 0 1; 0 0 0 -Bm/J];
B = [0 0; Kt*Ka/J 0; 0 0; 0 Kt*Ka/J];
C = [1 0 0 0; 0 0 1 0];
D = zeros(2,2);

sys_c = ss(A,B,C,D);

Ts =  0.01;

sys_d_fast = c2d(sys_c,Ts);

[Ad,Bd,Cd,Dd] = ssdata(sys_d_fast);

nx = 4;
nu = 2;
nr = 2;
ny = 2;

%% System fast
Ts_f =  0.01;

sys_d_fast = c2d(sys_c,Ts_f);

[Ad_f,Bd_f,Cd_f,Dd_f] = ssdata(sys_d_fast);

nx = 4;
nu = 2;
nr = 2;
ny = 2;


%% Augment States for reference tracking formulation 2
%x = [theta1; omega1; theta2; omega2; u1; u2; r1; r2]

Ax = [Ad                Bd          zeros(nx,nr)
      zeros(nu,nx)   eye(nu)        zeros(nu,nr)       
      zeros(nr,nx)   zeros(nr,nu)      eye(nr)];

Bx = [Bd;eye(nu);zeros(nr,nu)];

Ex = [Cd zeros(ny, nu) -eye(nr)];

Cx = eye(8);

Dx = zeros(8,2);

model = LTISystem('A', Ax , 'B', Bx, 'C', Cx, 'D', Dx, 'Ts', Ts);



%% Define constraints

% Control constraints
    model.u.min = [-1,-1]; %+-5V each motor to account for current limit (driver gain = 2);
    model.u.max= -model.u.min;

% State constraints
    model.x.min=[-deg2rad(40), -inf, -deg2rad(40), -inf, -5, -5, -inf, -inf];
    %model.x.min=[-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf];
    model.x.max=-model.x.min;
    
% Output constraints
    model.y.min = [-deg2rad(40), -inf, -deg2rad(40), -inf, -5, -5, -inf, -inf];
    %model.y.min=[-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf];
    model.y.max = -model.y.min;
    
%% Define penalties
Qy = 0.1;
Q = Ex'*Qy*Ex;
R = diag([0.1, 0.1]);
model.u.penalty = QuadFunction(R);
model.x.penalty = QuadFunction(Q);
%model.x.with('terminalPenalty');
%model.x.terminalPenalty = QuadFunction(Q/2);






%% Fast Controller M=2
N = 0.5/Ts;
N = 5;
ratio = Ts/0.01;
ctrl = MPCController(model,N);

%% Simulation setup
x0 = [0;0;0;0];
u0  = [0; 0];
r0 = [0; 0];

SimTime = floor(timeStop);

SimTime = 29;

Nsim = SimTime/Ts; %Simulation steps for entire IMU data ~20min


x = x0;
u = u0;
r =[0;0];
y =[0;0];

data.Time = ([1:Nsim]-1)*Ts;
data.X=zeros(nx,Nsim);
data.U=zeros(nu,Nsim);
data.R=zeros(nr,Nsim);
data.Y = zeros(2,Nsim);
data.Iter = zeros(1,Nsim);
data.J = zeros(1,Nsim);

%% Simulation

j = 0;
du_s = [0;0];
du_f = [0;0];
J = 0;

for i = 1:Nsim
    tic
    r(1) = theta1_d(i*ratio); %desired motor angles
    r(2) = theta2_d(i*ratio);
    
    [du,feasible,openloop] = ctrl.evaluate([x(:);u;r]);
    
    if isnan(du)
        du = [0;0];
        %du = -u;%-0.1*sign(u).*u;
    end
    
    loop_time = toc;
    delay_steps = round((loop_time - Ts)/Ts);
    %disp(loop_time);
    
    u = u + du;
    
    for j = 1:ratio
        x = Ad_f*x + Bd_f*u;
        y = Cd_f*x + Dd_f*u;
        xaug = [x; u; r];
        J = J + (xaug'*Q*xaug + du'*R*du);
    end
    
    
    %xaug = [x; u; r];
    %J = J + (xaug'*Q*xaug + du'*R*du)*ratio;
    
    data.X(:,i)=x;
    data.U(:,i)=u;
    data.R(:,i)=r;
    data.Y(:,i)=y;
    data.Iter(1,i) = loop_time;
    data.J(1,i) = J;
    
    
    
    if delay_steps>0 && i>10 && 0
        for j = 1:delay_steps
            x = Ad*x + Bd*u;
            y = Cd*x + Dd*u;
            data.X(:,i*j)=x;
            data.U(:,i*j)=u;
            data.R(:,i*j)=r;
            data.Y(:,i*j)=y;
            i = i+1;
            disp('hi');
        end
    end
    
    
end



%% Plot
figure(2)
plot(data.Time,data.X(1,:));
hold on
plot(data.Time,data.R(1,:));
legend('X1','R1');

figure(3)
plot(data.Time,data.U(1,:));
hold on
plot(data.Time,data.U(2,:));
legend('U1','U2');

figure(4)
plot(data.Time,data.U(1,:),'+');
hold on
%plot(data.Time_s,data.U_s(1,:),'o');


figure(5)
plot([1:Nsim],data.Iter(1,:).*1000,'+');
xlabel('Iteration');
ylabel('Computational time (ms)')
yline(5);

figure(6)
plot([1:Nsim],data.J(1,:),'+');
xlabel('Iteration');
ylabel('Cost Function')

%% Save values

time_01_N5 = sum(data.Iter(1,:));
cost_01_N5 = data.J(end);




%% Find forward kinematics using lookup table. Ie, Recover alpha and gamma

%{
load('lookup_angles.mat')
for i = 1:length(data.X(1,:))
    theta1 = data.X(1,i);
    theta2 = data.X(3,i);
    poss_ind = intersect(find(abs(lookup_angles(:,1)-theta1)<0.5),find(abs(lookup_angles(:,2)-theta2)<0.5));
    index = round(median(poss_ind));
    alpha_f = lookup_angles(index,3);
    gamma_f = lookup_angles(index,4);
end
%}
