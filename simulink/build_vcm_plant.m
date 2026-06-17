%% build_vcm_plant.m  - Model 1: VCM-유연기구 개루프 플랜트 (Simulink) + 개루프 해석
%  전기(R,L,역기전력) + 전자기(F=Kf*I) + 기계(m x'' + c x' + k x = F)
clear; clc; bdclose all;

%% 설계 파라미터 (권장 d_c=0.3 설계)
m=0.0147; k=6441; Kf=3.392; R=2.58; L=1.5e-3;
zeta=0.05; wn=sqrt(k/m); c=2*zeta*sqrt(k*m);
fprintf('wn=%.1f rad/s (f1=%.1f Hz), c=%.4f N s/m\n', wn, wn/2/pi, c);

mdl='vcm_plant'; try, close_system(mdl,0); catch, end; new_system(mdl);

add_block('simulink/Sources/Step',[mdl '/Vin'],'Position',[20 100 50 130],'Time','0','Before','0','After','1');
add_block('simulink/Math Operations/Sum',[mdl '/Sum'],'Position',[90 98 110 132],'Inputs','+-');
add_block('simulink/Continuous/Transfer Fcn',[mdl '/Elec'],'Position',[150 92 230 138],...
    'Numerator','[1]','Denominator',['[' num2str(L) ' ' num2str(R) ']']);
add_block('simulink/Math Operations/Gain',[mdl '/KfF'],'Position',[260 100 300 130],'Gain',num2str(Kf));
add_block('simulink/Continuous/Transfer Fcn',[mdl '/MechPos'],'Position',[340 92 440 138],...
    'Numerator','[1]','Denominator',['[' num2str(m) ' ' num2str(c) ' ' num2str(k) ']']);
add_block('simulink/Continuous/Transfer Fcn',[mdl '/MechVel'],'Position',[340 180 440 226],...
    'Numerator','[1 0]','Denominator',['[' num2str(m) ' ' num2str(c) ' ' num2str(k) ']']);
add_block('simulink/Math Operations/Gain',[mdl '/KfB'],'Position',[260 185 300 215],'Gain',num2str(Kf));
add_block('simulink/Sinks/To Workspace',[mdl '/Xout'],'Position',[480 100 540 130],...
    'VariableName','Xsig','SaveFormat','Timeseries');

add_line(mdl,'Vin/1','Sum/1');
add_line(mdl,'Sum/1','Elec/1');
add_line(mdl,'Elec/1','KfF/1');
add_line(mdl,'KfF/1','MechPos/1');
add_line(mdl,'MechPos/1','Xout/1');
add_line(mdl,'KfF/1','MechVel/1');     % 힘에서 분기 -> 속도
add_line(mdl,'MechVel/1','KfB/1');
add_line(mdl,'KfB/1','Sum/2');         % 역기전력 피드백

set_param(mdl,'Solver','ode45','StopTime','0.06','SaveOutput','on');
save_system(mdl,[pwd filesep 'vcm_plant.slx']);
fprintf('saved vcm_plant.slx\n');

%% 시뮬레이션 (1V 스텝)
so=sim(mdl); X=so.Xsig;
figure('Position',[100 100 480 320]);
plot(X.Time*1e3, X.Data*1e3,'LineWidth',1.6); grid on;
xlabel('time [ms]'); ylabel('displacement [mm]');
title('Open-loop step response (1 V input)');
exportgraphics(gcf,fullfile('..','figure','fig9_plant_step.png'),'Resolution',200);

%% 개루프 해석 (해석적 TF로 Bode) : X/V = Kf/[(ms^2+cs+k)(Ls+R)+Kf^2 s]
s=tf('s'); G = Kf/((m*s^2+c*s+k)*(L*s+R)+Kf^2*s);
figure('Position',[100 100 480 360]); bode(G); grid on;
title('Open-loop plant Bode (V \rightarrow x)');
exportgraphics(gcf,fullfile('..','figure','fig10_plant_bode.png'),'Resolution',200);
[gpeak,fpk]=getPeakGain(G); fprintf('공진 피크 @ %.1f Hz\n', fpk/2/pi);
fprintf('Model 1 완료: vcm_plant.slx + fig9/fig10\n');
