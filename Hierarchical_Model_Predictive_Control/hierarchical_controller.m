%clear all
%clc
%% Input IMU Data
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




%% Slow Controller M=1
Ts_slow =  0.1;

sys_d_slow = c2d(sys_c,Ts_slow);

[Ad_s,Bd_s,Cd_s,Dd_s] = ssdata(sys_d_slow);

Ax_s = [Ad_s                Bd_s          zeros(nx,nr)
      zeros(nu,nx)   eye(nu)        zeros(nu,nr)       
      zeros(nr,nx)   zeros(nr,nu)      eye(nr)];

Bx_s = [Bd_s;eye(nu);zeros(nr,nu)];

Ex_s = [Cd_s zeros(ny, nu) -eye(nr)];

Cx_s = eye(8);

Dx_s = zeros(8,2);

model_s = LTISystem('A', Ax_s , 'B', Bx_s, 'C', Cx_s, 'D', Dx_s, 'Ts', Ts_slow);

% Control constraints
    model_s.u.min= [-1,-1]; %+-5V each motor to account for current limit (driver gain = 2);
    model_s.u.max= -model.u.min;

% State constraints
    model_s.x.min=[-deg2rad(40), -inf, -deg2rad(40), -inf, -5, -5, -inf, -inf];
    %model.x.min=[-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf];
    model_s.x.max=-model.x.min;
    
% Output constraints
    model_s.y.min = [-deg2rad(40), -inf, -deg2rad(40), -inf, -5, -5, -inf, -inf];
    %model.y.min=[-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf];
    model_s.y.max = -model.y.min;

Qy = 1;
Qn = Ex'*Qy*Ex;
Rn = diag([0.1, 0.1]);
model_s.u.penalty = QuadFunction(Rn);
model_s.x.penalty = QuadFunction(Qn);
%model_s.x.with('terminalPenalty');
%model_s.x.terminalPenalty = QuadFunction(Q/2);

N_s = 0.5/Ts_slow;
ctrl_s = MPCController(model_s,N_s);

%% Fast Controller M=2
N = Ts_slow/Ts;
ctrl = MPCController(model,N);

%% Simulation setup
x0 = [0;0;0;0];
u0  = [0; 0];
r0 = [0; 0];

SimTime = floor(timeStop);

SimTime = 29;

Nsim = SimTime/Ts; %Simulation steps for entire IMU data ~20min

Nsim_s = SimTime/Ts_slow;

steps_ratio = Nsim/Nsim_s;

x = x0;
u = u0;
r =[0;0];
y =[0;0];
J = 0;

data.Time = ([1:Nsim]-1)*Ts;
data.Time_s = ([1:Nsim_s]-1)*Ts_slow;
data.X=zeros(nx,Nsim);
data.X_s=zeros(nx,Nsim_s);
data.U=zeros(nu,Nsim);
data.U_s=zeros(nu,Nsim_s);
data.R=zeros(nr,Nsim);
data.Y = zeros(2,Nsim);
data.Iter = zeros(1,Nsim);
data.J = zeros(1,Nsim);

%% Simulation

j = 0;
i = 1;
du_s = [0;0];
du_f = [0;0];
while j<Nsim_s
    i = 1;
    while i<=steps_ratio
    tic
        r(1) = theta1_d(i+j*steps_ratio); %desired motor angles
        r(2) = theta2_d(i+j*steps_ratio);
        r(1) = 0.8*sin((i+j*steps_ratio)/20); %desired motor angles
        r(2) = 0.8*cos((i+j*steps_ratio)/20);
        %try
            if i==1 %|| i==steps_ratio
                %r(1) = theta1_d(i+(j+1)*steps_ratio); %desired motor angles
                %r(2) = theta2_d(i+(j+1)*steps_ratio);
                %r(1) = sin((i+(j+1)*steps_ratio)/10); %desired motor angles
                %r(2) = cos((i+(j+1)*steps_ratio)/10);
                [du_s,feasible,openloop] = ctrl_s.evaluate([x(:);u;r]);
                if isnan(du_s)
                    du_s = [0;0];%-u;%-0.1*sign(u).*u;
                    %disp(j);
                end
                du_s = du_s; %change ratio
                %du = du_s;
                %u = u+du_s;
                du = du_s;
                data.U_s(:,j+1) = u+du_s;
                data.X_s(:,j+1) = Ad_s*x + Bd_s*(u+du_s); %predicted next state by slow controller
                %r(1) = theta1_d(i+j*steps_ratio); %desired motor angles
                %r(2) = theta2_d(i+j*steps_ratio);
                %r(1) = sin((i+j*steps_ratio)/10); %desired motor angles
                %r(2) = cos((i+j*steps_ratio)/10);
            else
                du_s = [0;0];
                [du,feasible,openloop] = ctrl.evaluate([x(:);u+du_s;r]);
            end
        %catch
            %du = 0;
        %end
        %[du,feasible,openloop] = ctrl.evaluate([x(:);u;r]);
        
        if isnan(du)
            du = [0;0];
            %du = -u;%-0.1*sign(u).*u;
            %disp(j);
        end

        loop_time = toc;  
        delay_steps = round((loop_time - Ts)/Ts);
        %disp(loop_time);
        
        u = u + du; % + du_s;
        x = Ad*x + Bd*u;
        y = Cd*x + Dd*u;
        xaug = [x; u; r];
        J = J + (xaug'*Q*xaug + du'*R*du);
        
        data.X(:,i+j*steps_ratio)=x;
        data.U(:,i+j*steps_ratio)=u;
        data.R(:,i+j*steps_ratio)=r;
        data.Y(:,i+j*steps_ratio)=y;
        data.J(:,i+j*steps_ratio)=J;
        data.Iter(:,i+j*steps_ratio) = loop_time;
        
        i = i+1;

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
    
    j = j+1;
end

%% Plot
figure(2)
stairs(data.Time,data.X(1,:),'LineWidth',1.5);
hold on
stairs(data.Time,data.R(1,:),'LineWidth',1.5);
hold on
stairs(data.Time,data.X(3,:),'LineWidth',1.5);
hold on
stairs(data.Time,data.R(2,:),'LineWidth',1.5);
legend('\theta_1','R1','\theta_2','R2');

figure(3)
plot(data.Time,data.U(1,:));
hold on
plot(data.Time,data.U(2,:));
legend('U1','U2');

figure(4)
stairs(data.Time,data.U(1,:),'LineWidth',1.5);
hold on
stairs(data.Time_s,data.U_s(1,:),'LineWidth',1.5);
axi=gca;
axi.FontSize = 13;
xlabel('Time(s)');
ylabel('u_1 (A)');
legend('Lower MPC', 'Upper MPC')
xlim([25,26.2])

figure(5)
stairs(data.Time,data.X(1,:),'LineWidth',1.5);
hold on
stairs(data.Time_s,data.X_s(1,:),'LineWidth',1.5);
hold on
plot(data.Time,data.R(1,:));
axi=gca;
axi.FontSize = 13;
xlabel('Time(s)');
ylabel('\theta_1 (rad)');
legend('Lower MPC', 'Upper MPC','Reference')
xlim([25,26.2])

figure(6)
plot([1:Nsim],data.Iter(1,:).*1000,'+');
xlabel('Iteration');
ylabel('Computational time (ms)')
yline(5);


%% Save values
time_hier = sum(data.Iter(1,:));
cost_hier = data.J(end);
save('hier_qyL100_qyU10_R001_sim29.mat','data');

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
