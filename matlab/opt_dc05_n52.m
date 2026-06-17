%% opt_dc05_n52.m
%  교수님 스펙 그대로(d_c>=0.5mm) 충족시키는 재최적화
%  변경점 (opt_design.m 대비):
%   - d_c = 0.5mm 고정 (스펙 준수)
%   - 영구자석 N52 (Br=1.43T): B_g surrogate를 Br비(1.43/1.17)로 선형 스케일
%   - 요크 FeCo(Hiperco): 중앙요크 포화한계 By1 1.8 -> 2.2T 완화
%   - 가동부 경량화(중공): 가동 stage/중간 stage 질량에 hollow 계수
%  설계변수 x=[t,L,b,r_m,h_m,h_y1,t_g, gamma] (mm), d_c는 고정
%  목적: gamma 최대화 (stroke>=gamma*1mm, f1>=gamma*100Hz). gamma>=1 이면 동시충족.
clear; clc; rng(7);

DC    = 0.5e-3;     % 코일 선경 고정 (스펙)
BR    = 1.43;       % N52 잔류자속 [T]
BR0   = 1.17;       % surrogate 적합시 자석 등급(N35)
HOLLOW= 0.4;        % 중공 경량화: solid 대비 잔존 질량비 (stage/중간stage)
BY1MAX= 2.2;        % FeCo 요크 포화 허용 [T]
AS    = 12e-3;      % 가동 stage 한 변 [m] (경량화: 15->12mm)

% 재료: 가동부+유연기구 (E[Pa], rho[kg/m3], Sy[Pa])
mats = { 'Al6061-T6', 71e9, 2770, 276e6;
         'Mg AZ31',   45e9, 1740, 200e6 };

% [ t    L    b    r_m  h_m  h_y1 t_g  gamma]
lb = [0.3, 20,  8,   2.0, 25,  4,   1.2, 0.2];
ub = [1.5, 45,  18,  6.5, 25,  16,  3.0, 3.0];
x0 = [0.6, 35,  15,  5.0, 25,  9,   2.0, 0.5];

opts = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
        'MaxIterations',300,'MaxFunctionEvaluations',3000, ...
        'StepTolerance',1e-11,'ConstraintTolerance',1e-9);

cfg.DC=DC; cfg.BR=BR; cfg.BR0=BR0; cfg.HOLLOW=HOLLOW; cfg.BY1MAX=BY1MAX; cfg.AS=AS;

fprintf('=== 교수님 스펙 준수 재최적화 (d_c=%.1fmm 고정, N52, FeCo요크, 중공%.0f%%) ===\n', ...
        DC*1e3, (1-HOLLOW)*100);
fprintf('%-12s%8s%11s%9s%10s%7s%8s%9s%9s\n', ...
    'Material','gamma','stroke[mm]','f1[Hz]','sig[MPa]','SF','Kf[N/A]','Bg[mT]','m_eff[g]');
fprintf('%s\n',repmat('-',1,84));
BX=cell(2,1); BR_=cell(2,1);
for k=1:2
    mat.name=mats{k,1}; mat.E=mats{k,2}; mat.rho=mats{k,3}; mat.Sy=mats{k,4};
    [bx,br]=solve_opt(lb,ub,x0,mat,cfg,opts,14); BX{k}=bx; BR_{k}=br;
    verdict=''; if bx(8)>=1.0, verdict=' <= 충족'; end
    fprintf('%-12s%8.3f%11.3f%9.1f%10.1f%7.2f%8.3f%9.1f%9.1f%s\n', ...
        mat.name, bx(8), br.x_max*1e3, br.f1, br.sigma*1e-6, br.SF, br.Kf, br.Bg*1e3, br.m_eff*1e3, verdict);
end

% --- 권장안 상세 (gamma 큰 쪽) ---
[~,kbest]=max([BX{1}(8),BX{2}(8)]); bx=BX{kbest}; br=BR_{kbest};
fprintf('\n--- 권장 설계 상세: %s ---\n', mats{kbest,1});
fprintf('유연기구: t=%.2f mm, L=%.1f mm, b=%.1f mm,  k_x=%.0f N/m\n', bx(1),bx(2),bx(3),br.k_x);
fprintf('VCM:      r_m=%.2f, h_m=%.1f, h_y1=%.2f, t_g=%.2f mm (d_c=%.1f), Ø=%.1f mm\n', ...
        bx(4),bx(5),bx(6),bx(7),DC*1e3, br.VCMd*1e3);
fprintf('자기:     Bg=%.0f mT, By1=%.2f T (<=%.1f), n=%.0f turns, m_coil=%.1f g\n', ...
        br.Bg*1e3, br.By1, BY1MAX, round(br.n), br.m_coil*1e3);
fprintf('성능:     f1=%.1f Hz, stroke@2A=%.3f mm, 1mm구동전류=%.2f A, SF=%.2f\n', ...
        br.f1, br.x_max*1e3, 1e-3/(br.Kf/br.k_x), br.SF);
fprintf('스펙판정: f1>100Hz [%s],  stroke>1mm [%s],  d_c>=0.5mm [OK],  Ø<=20mm [%s]\n', ...
        tf(br.f1>100), tf(br.x_max*1e3>1.0), tf(br.VCMd<=20e-3));

%% ===== 중공률 스캔 (Mg) : gamma=1 임계 경량화 찾기 =====
fprintf('\n=== 중공률 스캔 (Mg AZ31, a_s=%.0fmm) : gamma>=1 임계 ===\n', AS*1e3);
fprintf('%-12s%8s%11s%9s%9s\n','잔존질량비','gamma','stroke[mm]','f1[Hz]','판정');
fprintf('%s\n',repmat('-',1,50));
matM.name='Mg AZ31'; matM.E=45e9; matM.rho=1740; matM.Sy=200e6;
for hol=[0.5 0.4 0.3 0.2]
    cfg2=cfg; cfg2.HOLLOW=hol;
    [bx,br]=solve_opt(lb,ub,x0,matM,cfg2,opts,12);
    verdict='미충족'; if bx(8)>=1.0, verdict='충족 OK'; end
    fprintf('%-12.0f%8.3f%11.3f%9.1f%9s\n', hol*100, bx(8), br.x_max*1e3, br.f1, verdict);
end

%% ================== 함수 ==================
function s=tf(b), if b, s='OK'; else, s='X'; end, end

function [bestx,bestr]=solve_opt(lb,ub,x0,mat,cfg,opts,nstart)
    bestg=-Inf; bestx=x0;
    starts=[x0; lb+rand(nstart,8).*(ub-lb)];
    for s=1:size(starts,1)
        try
            [xo,fv,ef]=fmincon(@(x)-x(8), starts(s,:), [],[],[],[], lb,ub, @(x)nlcon(x,mat,cfg), opts);
            if ef>0
                c=nlcon(xo,mat,cfg);
                if all(c<=1e-6) && (-fv)>bestg, bestg=-fv; bestx=xo; end
            end
        catch
        end
    end
    bestr=perf(bestx,mat,cfg);
end

function r=perf(xv,mat,cfg)
    mu0=4*pi*1e-7;
    t=xv(1)*1e-3; L=xv(2)*1e-3; b=xv(3)*1e-3;
    r_m=xv(4)*1e-3; h_m=xv(5)*1e-3; h_y1=xv(6)*1e-3; t_g=xv(7)*1e-3; d_c=cfg.DC;
    t_y3=1e-3; h_y2=3e-3; c_c=0.3e-3; t_b=0.3e-3; r_p=1.15;
    n_sym=2; I_max=2; a_s=cfg.AS; d_b=40e-3; h_i=5e-3;
    w_c=t_g-2*c_c; h_c=h_y1;

    % B_g: FEMM-적합 surrogate(N35) -> N52로 Br비 선형 스케일
    Cs=[199.0491,92.0382,-30.2822,-82.8183,-0.3287,1.0510,8.0321,-3.1713,-7.3583,4.5247];
    rm=r_m*1e3; hy=h_y1*1e3; tg=t_g*1e3;
    Bg=(Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000 * (cfg.BR/cfg.BR0);
    A_g=2*pi*(r_m+t_g/2)*h_y1; Phi_g=Bg*A_g; By1=Phi_g/(pi*r_m^2);
    dcs=r_p*d_c;
    n_half=(h_c/dcs)*((2/sqrt(3))*(max(w_c,0)/dcs));
    L_m=2*pi*(r_m+c_c+t_b+w_c/2);
    n_eff=n_sym*(h_y1+0.08*h_m)/h_c*n_half;
    n_eff=min(n_eff, n_sym*n_half);
    Kf=n_eff*Bg*L_m;
    m_coil=8960*(n_sym*n_half)*(pi/4)*d_c^2*L_m;

    Ia=b*t^3/12; k_x=24*mat.E*Ia/L^3;
    m_ms =cfg.HOLLOW*mat.rho*a_s^2*b;          % 가동 stage(중공 경량화)
    m_int=cfg.HOLLOW*2*mat.rho*d_b*h_i*b;      % 중간 stage(중공)
    m_bld=8*mat.rho*L*t*b;                     % blade(전질량, 1/3 유효)
    m_eff=m_ms+m_coil+0.25*m_int+(1/3)*m_bld;
    f1=(1/(2*pi))*sqrt(k_x/m_eff);
    F_max=Kf*I_max; x_max=F_max/k_x;
    sigma=3*mat.E*t*(x_max/2)/L^2; SF=mat.Sy/sigma;
    VCMd=2*(r_m+t_g+t_y3); footY=a_s+2*(2*L+h_i);

    r.f1=f1; r.k_x=k_x; r.m_eff=m_eff; r.x_max=x_max; r.sigma=sigma; r.SF=SF;
    r.Kf=Kf; r.Bg=Bg; r.By1=By1; r.VCMd=VCMd; r.footY=footY; r.m_coil=m_coil;
    r.w_c=w_c; r.dcs=dcs; r.n=n_sym*n_half;
end

function [c,ceq]=nlcon(xv,mat,cfg)
    r=perf(xv,mat,cfg); g=xv(8); x_req=1e-3; f_req=100;
    c(1)=g*x_req - r.x_max;
    c(2)=g*f_req - r.f1;
    c(3)=r.sigma - 0.3*mat.Sy;
    c(4)=r.By1 - cfg.BY1MAX;
    c(5)=r.VCMd - 20e-3;
    c(6)=r.footY - 200e-3;
    c(7)=r.dcs - r.w_c;
    c(8)=r.Bg - 0.7;
    ceq=[];
end
