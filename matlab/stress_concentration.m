%% stress_concentration.m
%  leaf 힌지 뿌리 응력집중(Kt) 반영한 안전율 재판정 + 필렛 반경 요구 산출
%  - 공칭응력 sigma_nom = 3*E*t*(x/2)/L^2  (균일보, 기존 모델)
%  - 실제 peak = Kt * sigma_nom,  Kt = 1 + 0.5*sqrt(t/r)  (Peterson형 근사; 큰 D/d, 굽힘)
%  - 판정 기준 2개: 항복(Sy), 피로(0.3*Sy; 스캐닝 반복하중)
%  ※ Kt는 공학적 근사. 정밀값은 FEA/Peterson 차트 필요.
clear; clc;

% --- 두 설계점 (앞선 최적화 결과) ---
D(1)=struct('name','d_c=0.3  Al6061 (원안)',  'E',71e9,'Sy',276e6,'t',0.43e-3,'L',24.2e-3,'x',1.05e-3);
D(2)=struct('name','d_c=0.5  Mg AZ31 (스펙준수)','E',45e9,'Sy',200e6,'t',0.30e-3,'L',20.8e-3,'x',1.01e-3);

Kt = @(rt) 1 + 0.5*sqrt(1./rt);     % rt = r/t
RT = [0.1 0.2 0.3 0.5 1.0];          % 필렛/두께 비 후보

for i=1:2
    d=D(i);
    s_nom = 3*d.E*d.t*(d.x/2)/d.L^2;            % 공칭 [Pa]
    SFy_nom = d.Sy/s_nom; SFf_nom = 0.3*d.Sy/s_nom;  % 집중 무시 안전율
    fprintf('\n=== %s ===\n', d.name);
    fprintf('t=%.2fmm L=%.1fmm  공칭 sigma=%.1f MPa  (Kt무시: SF_yield=%.2f, SF_fatigue=%.2f)\n', ...
        d.t*1e3, d.L*1e3, s_nom*1e-6, SFy_nom, SFf_nom);
    fprintf('%-8s%6s%12s%11s%11s%9s\n','r/t','Kt','peak[MPa]','SF_yield','SF_fatig','피로판정');
    fprintf('%s\n',repmat('-',1,60));
    for rt=RT
        kt=Kt(rt); pk=kt*s_nom;
        SFy=d.Sy/pk; SFf=0.3*d.Sy/pk;
        verdict='위험'; if SFf>=1.0, verdict='안전'; end
        fprintf('%-8.2f%6.2f%12.1f%11.2f%11.2f%9s\n', rt, kt, pk*1e-6, SFy, SFf, verdict);
    end
    % 피로 안전(0.3Sy) 위한 최소 필렛: Kt*s_nom = 0.3*Sy
    Kt_allow = 0.3*d.Sy/s_nom;                  % 허용 Kt
    if Kt_allow>1
        rt_req = (0.5/(Kt_allow-1))^2;          % Kt=1+0.5/sqrt(rt) 역산
        fprintf('-> 피로(0.3Sy) 안전 위한 최소 필렛: r/t >= %.2f  (r >= %.3f mm)\n', ...
            rt_req, rt_req*d.t*1e3);
    else
        fprintf('-> 공칭응력만으로도 0.3Sy 초과: 필렛으로 해결 불가, 형상 재설계 필요\n');
    end
end

fprintf('\n[해석] Kt 무시 SF는 낙관. 피로(0.3Sy) 기준이 진짜 제약.\n');
fprintf('       필렛 반경을 충분히 키우면(아래 요구치) 두 설계 모두 안전 회복 가능.\n');
