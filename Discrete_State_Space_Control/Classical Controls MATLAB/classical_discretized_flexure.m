s = tf('s');
Gyy= -8.38e6*(s^2+21.12*s+2.28e4)*(s^2+11.27*s+1.58e6)*...
(s^2+3.25*s+4.72e6)*(s^2+20.72*s+6.71e6)*(s-7500)/...
    ((s^2+5.77*s+1.164e4)*(s^2+24.43*s+3.14e4)*...
    (s^2+17.50*s+1.54e6)*(s^2+6.99*s+4.49e6)*...
    (s^2+14.76*s+7.04e6)*(s+7500));
Cyy=260*(s+100)^2/(s*(s+2000)*(s+3000));
Lyy=Gyy*Cyy;
Tyy=Lyy/(1+Lyy);
Syy=1-Tyy;

Ts=0.0001;
Cyy=ss(Cyy);
Gyy=ss(Gyy);
Cyy=c2d(Cyy,Ts,'zoh');
Gyy=c2d(Gyy,Ts,'zoh');
CG_ssd = Cyy*Gyy;
T_ssd = feedback(CG_ssd,1);
S_ssd = 1-T_ssd;

%Controller designed using sisotool
z=tf('z');
Cyyz=0.14*(z-0.994)*(z-0.976)/((z-1)*(z+0.2)*(z-0.8));
CG_ssdz = Cyyz*Gyy;
T_ssdz = feedback(CG_ssdz,1);
S_ssdz = 1-T_ssdz;


figure
k=bodeplot(S_ssd);
title('Gaurav Discretized Controller Sensitivity');
setoptions(k,'FreqUnits','Hz','Xlim',[1 1000],'MagLowerLim',-40,'PhaseVisible','off');

figure
k=bodeplot(T_ssd);
title('Gaurav Discretized Controller Compl. Sensitivity');
setoptions(k,'FreqUnits','Hz','Xlim',[1 1000],'MagLowerLim',-40);

figure;
step(T_ssd);
title('Gaurav Discretized Controller Step Response');

figure
k=bodeplot(S_ssdz);
title('Designed Discretized Controller Sensitivity');
setoptions(k,'FreqUnits','Hz','Xlim',[1 1000],'MagLowerLim',-40,'PhaseVisible','off');

figure
k=bodeplot(T_ssd);
title('Designed Discretized Controller Compl. Sensitivity');
setoptions(k,'FreqUnits','Hz','Xlim',[1 1000],'MagLowerLim',-40);

figure;
step(T_ssdz);
title('Designed Discretized Controller Step Response');