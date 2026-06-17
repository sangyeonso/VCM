%% opt_kt_refinal.m
%  응력집중 K_t를 피로 제약에 반영한 최종 재최적화
%  제약 변경: 기존 sigma<=0.3Sy  ->  K_t*sigma_nom <= 0.3Sy  (피로, 반복하중)
%  K_t=1.7 (필렛 r/t≈0.5, Peterson형 근사). 그 외 물리/VCM 모델은 동일.
%  설계변수 x=[t,L,b,r_m,h_m,h_y1,t_g, gamma] (mm), d_c는 cfg로 고정
clear; clc; rng(7);

opts = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
        'MaxIterations',300,'MaxFunctionEvaluations',3000, ...
        'StepTolerance',1e-11,'ConstraintTolerance',1e-9);
lb = [0.3, 20,  8,   2.0, 25,  4,   1.0, 0.2];
ub = [1.5, 45,  18,  6.5, 25,  16,  3.0, 3.0];
x0 = [0.6, 35,  15,  5.0, 25,  9,   2.0, 0.5];

% 공통 cfg 기본값
base.BR0=1.17; base.KT=1.7;

% --- 트랙 정의 ---
%  A: d_c=0.3 원안 (Al, N35, solid, AISI요크, stage15mm)
%  B: d_c=0.5 스펙준수 (Mg, N52, 70%중공, FeCo요크, stage12mm)
TA=base; TA.DC=0.3e-3; TA.BR=1.17; TA.HOLLOW=1.0; TA.BY1MAX=1.8; TA.AS=15e-3;  % solid
matA.name='Al6061'; matA.E=71e9; matA.rho=2770; matA.Sy=276e6;
TB=base; TB.DC=0.5e-3; TB.BR=1.43; TB.HOLLOW=1.0; TB.BY1MAX=2.2; TB.AS=12e-3;
matB.name='Mg AZ31'; matB.E=45e9; matB.rho=1740; matB.Sy=200e6;

fprintf('=== K_t=%.1f (필렛 r/t≈0.5) 반영 재최적화 : 피로 K_t*sigma<=0.3Sy ===\n\n', base.KT);
runtrack('A: d_c=0.3 (Al, N35, 중공없음 solid)',  TA, matA, lb,ub,x0, opts);
runtrack('B: d_c=0.5 (Mg, N52, 중공없음 solid)',  TB, matB, lb,ub,x0, opts);

% --- K_t 민감도 (트랙 B) : 필렛 키우면 회복되나 ---
fprintf('\n=== K_t 민감도 (트랙 B, d_c=0.5) ===\n');
fprintf('%-22s%8s%11s%9s\n','K_t (필렛)','gamma','stroke[mm]','f1[Hz]');
fprintf('%s\n',repmat('-',1,52));
for kt=[2.1 1.7 1.5 1.3 1.0]
    cfg=TB; cfg.KT=kt;
    [bx,br]=solve_opt(lb,ub,x0,matB,cfg,opts,12);
    lbl=sprintf('%.1f (r/t=%.2f)',kt, (0.5/(kt-1))^2);
    if kt==1.0, lbl='1.0 (집중무시)'; end
    fprintf('%-22s%8.3f%11.3f%9.1f\n', lbl, bx(8), br.x_max*1e3, br.f1);
end
fprintf('\n주: gamma>=1 이면 stroke>=1mm & f1>=100Hz 동시충족(피로+K_t 안전 포함).\n');

%% ================== 함수 ==================
function runtrack(label,cfg,mat,lb,ub,x0,opts)
    [bx,br]=solve_opt(lb,ub,x0,mat,cfg,opts,16);
    sf_fat = 0.3*mat.Sy/(cfg.KT*br.sigma);   % 피로 안전율(K_t 포함)
    bind=''; if abs(cfg.KT*br.sigma-0.3*mat.Sy)<1e6, bind='[응력제약 binding]'; end
    fprintf('--- %s ---\n', label);
    fprintf('  gamma=%.3f  stroke@2A=%.3f mm  f1=%.1f Hz  %s\n', bx(8), br.x_max*1e3, br.f1, ...
        ternary(bx(8)>=1.0,'<= 스펙충족','<= 미충족'));
    fprintf('  공칭sig=%.1f MPa, peak(K_t)=%.1f MPa, 피로SF=%.2f %s\n', ...
        br.sigma*1e-6, cfg.KT*br.sigma*1e-6, sf_fat, bind);
    fprintf('  설계: t=%.2f L=%.1f b=%.1f mm | r_m=%.2f h_y1=%.2f t_g=%.2f | Bg=%.0fmT Kf=%.2f m_eff=%.1fg\n\n', ...
        bx(1),bx(2),bx(3),bx(4),bx(6),bx(7),br.Bg*1e3,br.Kf,br.m_eff*1e3);
end

function s=ternary(c,a,b), if c, s=a; else, s=b; end, end

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
    m_ms =cfg.HOLLOW*mat.rho*a_s^2*b;
    m_int=cfg.HOLLOW*2*mat.rho*d_b*h_i*b;
    m_bld=8*mat.rho*L*t*b;
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
    c(3)=cfg.KT*r.sigma - 0.3*mat.Sy;   % 피로 + 응력집중
    c(4)=r.By1 - cfg.BY1MAX;
    c(5)=r.VCMd - 20e-3;
    c(6)=r.footY - 200e-3;
    c(7)=r.dcs - r.w_c;
    c(8)=r.Bg - 0.7;
    ceq=[];
end
