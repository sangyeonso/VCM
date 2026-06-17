%% build_vcm_simscape.m - Model 3: VCM-유연기구 Simscape 물리네트워크 (완전 자동연결+검증)
%  결합을 Simulink 도메인으로 우회(PS->Simulink->xKf->Simulink->PS, 단위 N/V 지정)해 단위정합.
clear; clc; bdclose all;
m=0.0147; k=6441; Kf=3.392; R=2.58; L=1.5e-3; zeta=0.05; c=2*zeta*sqrt(k*m);
NL=char(10); P=@(s) strrep(s,'|',NL);
for lib={'fl_lib','nesl_utility','simulink'}, try, load_system(lib{1}); catch, end; end
mdl='vcm_simscape'; try, close_system(mdl,0); catch, end; new_system(mdl);

BLK={
 'nesl_utility/Solver Configuration','Solver',[60 300 110 330],{};
 'nesl_utility/Solver Configuration','Solver2',[700 450 750 480],{};
 'fl_lib/Electrical/Electrical Elements/Electrical Reference','GND_e',[60 360 95 395],{};
 'fl_lib/Electrical/Electrical Sources/Controlled Voltage|Source','Vsrc',[170 180 210 240],{};
 'fl_lib/Electrical/Electrical Elements/Resistor','Rcoil',[250 175 305 205],{'R',num2str(R)};
 'fl_lib/Electrical/Electrical Elements/Inductor','Lcoil',[330 175 385 205],{'l',num2str(L)};
 'fl_lib/Electrical/Electrical Sensors/Current Sensor','Isens',[410 180 460 230],{};
 'fl_lib/Electrical/Electrical Sources/Controlled Voltage|Source','Vbemf',[170 290 210 350],{};
 % --- 전류->힘 결합 (Simulink 우회) ---
 'nesl_utility/PS-Simulink Converter','P2S_I',[490 120 525 150],{'Unit','A'};
 'simulink/Math Operations/Gain','GainF',[545 120 585 150],{'Gain',num2str(Kf)};
 'nesl_utility/Simulink-PS Converter','S2P_F',[605 120 640 150],{'Unit','N'};
 'fl_lib/Mechanical/Mechanical Sources/Ideal Force Source','Fsrc',[660 110 700 170],{};
 % --- 속도->역기전력 결합 ---
 'nesl_utility/PS-Simulink Converter','P2S_v',[300 410 335 440],{'Unit','m/s'};
 'simulink/Math Operations/Gain','GainB',[245 410 285 440],{'Gain',num2str(-Kf)};
 'nesl_utility/Simulink-PS Converter','S2P_V',[180 410 215 440],{'Unit','V'};
 % --- 기계 ---
 'fl_lib/Mechanical/Translational|Elements/Mass','Mass',[750 200 800 240],{'mass',num2str(m)};
 'fl_lib/Mechanical/Translational|Elements/Translational Spring','Spring',[750 290 800 330],{'spr_rate',num2str(k)};
 'fl_lib/Mechanical/Translational|Elements/Translational Damper','Damper',[830 290 880 330],{'D',num2str(c)};
 'fl_lib/Mechanical/Translational|Elements/Mechanical|Translational|Reference','GND_m',[770 390 805 425],{};
 'fl_lib/Mechanical/Mechanical Sensors/Ideal Translational|Motion Sensor','Msens',[830 180 880 240],{};
 % --- I/O ---
 'simulink/Sources/Step','Vin',[40 130 70 160],{'Time','0.005','Before','0','After','1'};
 'nesl_utility/Simulink-PS Converter','S2PS',[100 130 135 160],{'Unit','V'};
 'nesl_utility/PS-Simulink Converter','PS2S',[930 190 965 220],{'Unit','m'};
 'simulink/Sinks/To Workspace','Xo',[1000 190 1050 220],{'VariableName','Xss','SaveFormat','Timeseries'};
};
for i=1:size(BLK,1), add_block(P(BLK{i,1}),[mdl '/' BLK{i,2}],'Position',BLK{i,3},BLK{i,4}{:}); end

LN={ {'Vin/1','S2PS/1'},{'S2PS/RConn1','Vsrc/RConn1'}, ...
 {'Vsrc/LConn1','Rcoil/LConn1'},{'Rcoil/RConn1','Lcoil/LConn1'}, ...
 {'Lcoil/RConn1','Isens/LConn1'},{'Isens/RConn2','Vbemf/LConn1'}, ...
 {'Vbemf/RConn2','GND_e/LConn1'},{'Vsrc/RConn2','GND_e/LConn1'},{'Solver/RConn1','Vsrc/LConn1'}, ...
 {'Isens/RConn1','P2S_I/LConn1'},{'P2S_I/1','GainF/1'},{'GainF/1','S2P_F/1'},{'S2P_F/RConn1','Fsrc/RConn1'}, ...
 {'Msens/RConn2','P2S_v/LConn1'},{'P2S_v/1','GainB/1'},{'GainB/1','S2P_V/1'},{'S2P_V/RConn1','Vbemf/RConn1'}, ...
 {'Fsrc/RConn2','Mass/LConn1'},{'Fsrc/LConn1','GND_m/LConn1'}, ...
 {'Mass/LConn1','Spring/LConn1'},{'Spring/RConn1','GND_m/LConn1'}, ...
 {'Mass/LConn1','Damper/LConn1'},{'Damper/RConn1','GND_m/LConn1'}, ...
 {'Mass/LConn1','Msens/LConn1'},{'Msens/RConn1','GND_m/LConn1'}, ...
 {'Solver2/RConn1','GND_m/LConn1'}, ...
 {'Msens/RConn3','PS2S/LConn1'},{'PS2S/1','Xo/1'} };
lf={}; for i=1:numel(LN), try, add_line(mdl,LN{i}{1},LN{i}{2},'autorouting','on'); catch, lf{end+1}=sprintf('%s->%s',LN{i}{1},LN{i}{2}); end; end %#ok
set_param(mdl,'Solver','ode23t','StopTime','0.08');
save_system(mdl,[pwd filesep 'vcm_simscape.slx']);
fprintf('연결 %d/%d 성공.\n',numel(LN)-numel(lf),numel(LN)); for i=1:numel(lf), fprintf(' 실패 %s\n',lf{i}); end

try
  so=sim(mdl); X=so.Xss; xss=X.Data(end)*1e3; xpk=max(X.Data)*1e3;
  fprintf('Simscape 시뮬 성공: 정상상태 x=%.3f mm (이론 %.3f mm), 피크 %.3f mm\n',...
    xss, (Kf/R)/k*1e3, xpk);
catch e, fprintf('시뮬 오류: %s\n', regexprep(getReport(e,'basic'),'\s+',' ')); end
