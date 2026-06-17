%% make_dynamics_figs.m - 동특성 figure 최신화 (최종 설계값, 전압모드)
%  전압모드: 전기(R,L,역기전력) 포함 -> 전기역학 감쇠로 ζ_eff≈0.99 (잘 감쇠)
%  fig11 폐루프 스텝(전류 모니터), fig12 사인 스캐닝.
clear; clc; close all;
Kf=1.04; Imax=2; f1=92.6;
k = Kf*Imax/0.93e-3;                 % 0.93mm@2A 정합
m = k/(2*pi*f1)^2;                   % f1=92.6Hz
R=0.15; L=0.1e-3; zeta=0.05; c=2*zeta*sqrt(k*m);
A=[0 1 0; -k/m -c/m Kf/m; 0 -Kf/L -R/L]; Bv=[0;0;1/L];
Gv=ss(A,Bv,[1 0 0],0);               % V -> x
Gi=ss(A,Bv,[0 0 1],0);               % V -> I
ze=(R*c+Kf^2)/(R*m)/(2*sqrt(k/m));
fprintf('전압모드 플랜트: ζ_eff=%.2f (전기역학 감쇠 K_f^2/R=%.1f)\n', ze, Kf^2/R);

C=pidtune(Gv,'PIDF');  CL=feedback(C,Gv);  T=feedback(C*Gv,1);
fprintf('폐루프 안정? %d\n', isstable(T));

%% fig11: 0.5mm 스텝 (변위 + 전류)
t=0:2e-5:0.05; ref=0.5e-3;
V=lsim(CL,ref*ones(size(t)),t);  y=lsim(Gv,V,t);  I=lsim(Gi,V,t);
figure('Position',[100 100 480 330]);
yyaxis left;  plot(t*1e3,y*1e3,'-','LineWidth',1.8); ylabel('변위 [mm]'); ylim([0 0.6]);
yline(0.5,'--','Color',[0.6 0.6 0.6]);
yyaxis right; plot(t*1e3,I,'-','LineWidth',1.2); ylabel('전류 [A]'); ylim([0 2.2]);
yline(Imax,':','Color',[0.8 0.2 0.2]);
xlabel('time [ms]'); grid on; title(sprintf('폐루프 0.5 mm 스텝 (전압모드, 최대전류 %.1f A)',max(abs(I))));
exportgraphics(gcf,fullfile('..','figure','fig11_step.png'),'Resolution',200);

%% fig12: 10Hz 사인 스캐닝 ±0.3mm
t=0:2e-5:0.3; fscan=10; ref=0.3e-3*sin(2*pi*fscan*t);
V=lsim(CL,ref,t); y=lsim(Gv,V,t);
err=ref(:)-y(:); rms_um=sqrt(mean(err.^2))*1e6;
figure('Position',[100 100 480 330]);
plot(t*1e3,ref*1e3,'--','LineWidth',1.3,'Color',[0.6 0.6 0.6]); hold on;
plot(t*1e3,y*1e3,'-','LineWidth',1.6,'Color',[0.1 0.3 0.9]);
xlabel('time [ms]'); ylabel('변위 [mm]'); grid on;
legend({'기준 (10Hz, \pm0.3mm)','추종'},'Location','northeast','FontSize',8);
title(sprintf('10 Hz 스캐닝 추종 (RMS 오차 %.0f \\mum)', rms_um));
exportgraphics(gcf,fullfile('..','figure','fig12_scan.png'),'Resolution',200);
fprintf('저장: fig11_step.png, fig12_scan.png (전압모드, 스텝 최대전류 %.2fA, 스캔 RMS=%.0fum)\n', max(abs(lsim(Gi,lsim(CL,ref,t),t))), rms_um);
