%% build_vcm_control.m - Model 2: 폐루프 위치제어 (PID + 2A 전류포화)
clear; clc; bdclose all;
m=0.0147; k=6441; Kf=3.392; R=2.58; L=1.5e-3; zeta=0.05; c=2*zeta*sqrt(k*m);

% 플랜트 V->x 로 PID 자동튜닝 (대역폭 ~ 60Hz, 공진 105Hz 아래)
s=tf('s'); G=Kf/((m*s^2+c*s+k)*(L*s+R)+Kf^2*s);
C=pidtune(G,'PIDF',2*pi*50);
if C.Tf>0, Nval=1/C.Tf; else, Nval=100; end   % Kd=0이면 N 임의(유한)
fprintf('PID: Kp=%.3g Ki=%.3g Kd=%.3g N=%.3g\n',C.Kp,C.Ki,C.Kd,Nval);

mdl='vcm_control'; try, close_system(mdl,0); catch, end; new_system(mdl);
add_block('simulink/Sources/Step',[mdl '/Ref'],'Position',[20 100 50 130],'Time','0.005','Before','0','After','0.001');
add_block('simulink/Math Operations/Sum',[mdl '/Err'],'Position',[80 98 100 132],'Inputs','+-');
add_block('simulink/Continuous/PID Controller',[mdl '/PID'],'Position',[130 95 200 135],...
    'P',num2str(C.Kp),'I',num2str(C.Ki),'D',num2str(C.Kd),'N',num2str(Nval));
add_block('simulink/Math Operations/Sum',[mdl '/Sb'],'Position',[230 98 250 132],'Inputs','+-');
add_block('simulink/Continuous/Transfer Fcn',[mdl '/Elec'],'Position',[270 92 350 138],...
    'Numerator','[1]','Denominator',['[' num2str(L) ' ' num2str(R) ']']);
add_block('simulink/Discontinuities/Saturation',[mdl '/Isat'],'Position',[370 100 400 130],...
    'UpperLimit','2','LowerLimit','-2');   % 2A 전류 한계
add_block('simulink/Math Operations/Gain',[mdl '/KfF'],'Position',[420 100 455 130],'Gain',num2str(Kf));
add_block('simulink/Continuous/Transfer Fcn',[mdl '/MechPos'],'Position',[480 92 580 138],...
    'Numerator','[1]','Denominator',['[' num2str(m) ' ' num2str(c) ' ' num2str(k) ']']);
add_block('simulink/Continuous/Transfer Fcn',[mdl '/MechVel'],'Position',[480 185 580 231],...
    'Numerator','[1 0]','Denominator',['[' num2str(m) ' ' num2str(c) ' ' num2str(k) ']']);
add_block('simulink/Math Operations/Gain',[mdl '/KfB'],'Position',[420 190 455 220],'Gain',num2str(Kf));
add_block('simulink/Sinks/To Workspace',[mdl '/Xo'],'Position',[620 100 670 130],'VariableName','Xsig','SaveFormat','Timeseries');
add_block('simulink/Sinks/To Workspace',[mdl '/Ro'],'Position',[620 40 670 70],'VariableName','Rsig','SaveFormat','Timeseries');
add_block('simulink/Sinks/To Workspace',[mdl '/Io'],'Position',[420 40 470 70],'VariableName','Isig','SaveFormat','Timeseries');

add_line(mdl,'Ref/1','Err/1'); add_line(mdl,'Ref/1','Ro/1');
add_line(mdl,'Err/1','PID/1'); add_line(mdl,'PID/1','Sb/1');
add_line(mdl,'Sb/1','Elec/1'); add_line(mdl,'Elec/1','Isat/1');
add_line(mdl,'Isat/1','KfF/1'); add_line(mdl,'Isat/1','Io/1');
add_line(mdl,'KfF/1','MechPos/1'); add_line(mdl,'MechPos/1','Xo/1');
add_line(mdl,'MechPos/1','Err/2');     % 위치 피드백
add_line(mdl,'KfF/1','MechVel/1'); add_line(mdl,'MechVel/1','KfB/1'); add_line(mdl,'KfB/1','Sb/2');

set_param(mdl,'Solver','ode23t','StopTime','0.08');
save_system(mdl,[pwd filesep 'vcm_control.slx']);

%% (a) 1mm 스텝응답
so=sim(mdl); X=so.Xsig; Rr=so.Rsig; I=so.Isig;
figure('Position',[100 100 480 330]);
plot(Rr.Time*1e3,Rr.Data*1e3,'--k','LineWidth',1.3); hold on;
plot(X.Time*1e3,X.Data*1e3,'b','LineWidth',1.6); grid on;
xlabel('time [ms]'); ylabel('position [mm]'); legend('reference','output','Location','southeast');
title('Closed-loop 1 mm step response');
info=stepinfo(X.Data,X.Time,0.001);
fprintf('정착시간 %.1f ms, 오버슈트 %.1f%%, 최대전류 %.2f A\n',info.SettlingTime*1e3,info.Overshoot,max(abs(I.Data)));
exportgraphics(gcf,fullfile('..','figure','fig11_step.png'),'Resolution',200);

%% (b) 10Hz 정현파 스캐닝 추종 (Ref -> Sine로 교체)
delete_line(mdl,'Ref/1','Err/1'); delete_line(mdl,'Ref/1','Ro/1');
delete_block([mdl '/Ref']);
add_block('simulink/Sources/Sine Wave',[mdl '/Ref'],'Position',[20 100 50 130],...
    'Amplitude','0.0005','Frequency','2*pi*10');
add_line(mdl,'Ref/1','Err/1'); add_line(mdl,'Ref/1','Ro/1');
set_param(mdl,'StopTime','0.3');
so=sim(mdl); X=so.Xsig; Rr=so.Rsig;
err=(interp1(X.Time,X.Data,Rr.Time)-Rr.Data); rms_err=sqrt(mean(err(Rr.Time>0.1).^2));
figure('Position',[100 100 480 330]);
plot(Rr.Time*1e3,Rr.Data*1e3,'--k','LineWidth',1.3); hold on;
plot(X.Time*1e3,X.Data*1e3,'b','LineWidth',1.4); grid on;
xlabel('time [ms]'); ylabel('position [mm]'); legend('reference','output','Location','northeast');
title(sprintf('10 Hz scanning tracking (RMS err %.0f \\mum)',rms_err*1e6));
exportgraphics(gcf,fullfile('..','figure','fig12_scan.png'),'Resolution',200);
save_system(mdl,[pwd filesep 'vcm_control.slx']);
fprintf('Model 2 완료: vcm_control.slx + fig11/fig12. 추종 RMS=%.0f um\n',rms_err*1e6);
