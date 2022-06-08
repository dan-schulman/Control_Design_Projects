%% Qy_L

QyL100 = load('hier_qyL100_qyU10_R01_sim29.mat','data');
QyL10 = load('hier_qyL10_qyU10_R01_sim29.mat','data');
QyL01 = load('hier_qyL01_qyU10_R01_sim29.mat','data');

figure(1)
stairs(QyL100.data.Time,QyL100.data.U(1,:),'LineWidth',1.5);
hold on
stairs(QyL10.data.Time,QyL10.data.U(1,:),'LineWidth',1.5);
hold on
stairs(QyL01.data.Time,QyL01.data.U(1,:),'LineWidth',1.5);
axi=gca;
axi.FontSize = 13;
xlabel('Time(s)');
ylabel('u_1 (A)');
legend('Qe^L = 10', 'Qe^L = 1', 'Qe^L = 0.1');
xlim([25.2,25.4])

figure(2)
stairs(QyL100.data.Time,QyL100.data.X(1,:),'LineWidth',1.5);
hold on
stairs(QyL10.data.Time,QyL10.data.X(1,:),'LineWidth',1.5);
hold on
stairs(QyL01.data.Time,QyL01.data.X(1,:),'LineWidth',1.5);
hold on
plot(QyL100.data.Time,QyL100.data.R(1,:));
axi=gca;
axi.FontSize = 13;
xlabel('Time(s)');
ylabel('\theta_1 (rad)');
legend('Qe^L = 10', 'Qe^L = 1', 'Qe^L = 0.1','Reference');
xlim([25.2,25.4])

%% Qy_U
QyU100 = load('hier_qyL10_qyU100_R01_sim29.mat','data');
QyU10 = load('hier_qyL10_qyU10_R01_sim29.mat','data');
QyU01 = load('hier_qyL10_qyU01_R01_sim29.mat','data');

figure(3)
stairs(QyU100.data.Time,QyU100.data.U(1,:),'LineWidth',1.5);
hold on
stairs(QyU10.data.Time,QyU10.data.U(1,:),'LineWidth',1.5);
hold on
stairs(QyU01.data.Time,QyU01.data.U(1,:),'LineWidth',1.5);
axi=gca;
axi.FontSize = 13;
xlabel('Time(s)');
ylabel('u_1 (A)');
legend('Qe^U = 10', 'Qe^U = 1', 'Qe^U = 0.1');
xlim([25.2,25.4])

figure(4)
stairs(QyU100.data.Time,QyU100.data.X(1,:)-0.005,'LineWidth',1.5);
hold on
stairs(QyU10.data.Time,QyU10.data.X(1,:),'LineWidth',1.5);
hold on
stairs(QyU01.data.Time,QyU01.data.X(1,:),'LineWidth',1.5);
hold on
plot(QyU100.data.Time,QyU100.data.R(1,:));
axi=gca;
axi.FontSize = 13;
xlabel('Time(s)');
ylabel('\theta_1 (rad)');
legend('Qe^U = 10', 'Qe^U = 1', 'Qe^U = 0.1','Reference');
xlim([25.2,25.4])


%% R

R01 = load('hier_qyL10_qyU10_R01_sim29.mat','data');
R001 = load('hier_qyL10_qyU10_R001_sim29.mat','data');
R0001 = load('hier_qyL10_qyU10_R0001_sim29.mat','data');

figure(5)
stairs(R01.data.Time,R01.data.U(1,:),'LineWidth',1.5);
hold on
stairs(R001.data.Time,R001.data.U(1,:),'LineWidth',1.5);
hold on
stairs(R0001.data.Time,R0001.data.U(1,:),'LineWidth',1.5);
axi=gca;
axi.FontSize = 13;
xlabel('Time(s)');
ylabel('u_1 (A)');
legend('R = 0.1', 'R = 0.01', 'R = 0.001');
xlim([25.2,25.4])

figure(6)
stairs(R01.data.Time,R01.data.X(1,:),'LineWidth',1.5);
hold on
stairs(R001.data.Time,R001.data.X(1,:),'LineWidth',1.5);
hold on
stairs(R0001.data.Time,R0001.data.X(1,:),'LineWidth',1.5);
hold on
plot(R01.data.Time,R01.data.R(1,:));
axi=gca;
axi.FontSize = 13;
xlabel('Time(s)');
ylabel('\theta_1 (rad)');
legend('R = 0.1', 'R = 0.01', 'R = 0.001','Reference');
xlim([25.2,25.4])