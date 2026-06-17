%% opt_hm_free.m
%  자석 길이 h_m을 envelope 한도(최대 100mm)까지 풀어 d_c=0.5 회복되는지 검증
%  - surrogate는 h_m 미포함 -> 해석적 퍼미언스로 B_g(h_m)/B_g(25) 비율 R만 곱함
%    R = (P_m(25)+P_l+P_main)/(P_m(h_m)+P_l+P_main),  P_m∝1/h_m
%  - 나머지는 no-hollow Track B (Mg, N52, FeCo, solid, a_s=12, K_t=1.7) 동일
clear; clc; rng(7);

cfg.DC=0.5e-3; cfg.BR=1.43; cfg.BR0=1.17; cfg.HOLLOW=1.0; cfg.BY1MAX=2.2;
cfg.AS=12e-3; cfg.KT=1.7;
mat.name='Mg AZ31'; mat.E=45e9; mat.rho=1740; mat.Sy=200e6;

opts = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
        'MaxIterations',300,'MaxFunctionEvaluations',3000, ...
        'StepTolerance',1e-11,'ConstraintTolerance',1e-9);
% [ t    L    b    r_m  h_m   h_y1 t_g  gamma]   <-- h_m 자유 [25,100]
lb = [0.3, 20,  8,   2.0, 25,   4,   1.0, 0.2];
ub = [1.5, 45,  18,  6.5, 100,  16,  3.0, 3.0];
x0 = [0.4, 26,  8,   6.5, 60,   4.5, 2.0, 0.9];

%% (1) h_m 단독 민감도 : Track B 최적형상 고정, h_m만 스윕
fprintf('=== (1) h_m 단독 효과 (형상: r_m=6.5,h_y1=4.5,t_g=2.0 고정) ===\n');
fprintf('%-9s%9s%9s%11s%9s\n','h_m[mm]','R(배율)','Bg[mT]','stroke[mm]','f1[Hz]');
fprintf('%s\n',repmat('-',1,48));
xfix=[0.4,26,8,6.5,25,4.5,2.0,0.9];
for hm=[25 40 60 80 100]
    xx=xfix; xx(5)=hm; r=perf(xx,mat,cfg);
    fprintf('%-9d%9.3f%9.0f%11.3f%9.1f\n', hm, r.R, r.Bg*1e3, r.x_max*1e3, r.f1);
end

%% (2) h_m 자유 재최적화
fprintf('\n=== (2) h_m 자유 재최적화 (d_c=0.5, Mg, N52, solid) ===\n');
[bx,br]=solve_opt(lb,ub,x0,mat,cfg,opts,18);
fprintf('gamma=%.3f  stroke@2A=%.3f mm  f1=%.1f Hz  %s\n', bx(8), br.x_max*1e3, br.f1, ...
    tern(bx(8)>=1.0,'<= 스펙충족','<= 미충족'));
fprintf('선택 h_m=%.1f mm (R=%.3f), Bg=%.0f mT, By1=%.2f T, Kf=%.2f, m_eff=%.1f g\n', ...
    bx(5), br.R, br.Bg*1e3, br.By1, br.Kf, br.m_eff*1e3);
fprintf('blade: t=%.2f L=%.1f b=%.1f mm | r_m=%.2f h_y1=%.2f t_g=%.2f mm\n', ...
    bx(1),bx(2),bx(3),bx(4),bx(6),bx(7));
fprintf('피로: 공칭sig=%.1f, peak(Kt)=%.1f MPa, SF=%.2f | VCM축길이≈h_m+요크=%.0f mm (<200 footprint)\n', ...
    br.sigma*1e-6, cfg.KT*br.sigma*1e-6, 0.3*mat.Sy/(cfg.KT*br.sigma), bx(5)+3);

%% ===== 함수 =====
function s=tern(c,a,b), if c, s=a; else, s=b; end, end

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
    Hc=867e3; mu_s=1400; n_sym=2; I_max=2; a_s=cfg.AS; d_b=40e-3; h_i=5e-3;
    w_c=t_g-2*c_c; h_c=h_y1;
    mu_r=cfg.BR/(mu0*Hc);

    % --- 해석적 퍼미언스로 h_m 효과 R 산출 (surrogate 보정용) ---
    P_y1=pi*mu0*mu_s*r_m^2/h_y1;
    P_y2=2*pi*mu0*mu_s*h_y2/(log(r_m+t_g+t_y3/2)-log(r_m/2));
    P_y3=pi*mu0*mu_s*((r_m+t_g+t_y3)^2-(r_m+t_g)^2)/(h_y1+h_m);
    P_g =2*pi*mu0*h_y1/(log(r_m+t_g)-log(r_m));
    P_l =0.52*(2*pi*mu0*r_m);
    P_main=1/(1/P_g+1/P_y1+1/P_y3+1/P_y2);
    P_m   =pi*mu0*mu_r*r_m^2/h_m;
    P_m25 =pi*mu0*mu_r*r_m^2/25e-3;
    R = (P_m25+P_l+P_main)/(P_m+P_l+P_main);   % h_m=25 대비 B_g 배율

    % --- surrogate(N35,h_m=25) -> N52 스케일 -> h_m 보정 R ---
    Cs=[199.0491,92.0382,-30.2822,-82.8183,-0.3287,1.0510,8.0321,-3.1713,-7.3583,4.5247];
    rm=r_m*1e3; hy=h_y1*1e3; tg=t_g*1e3;
    Bg=(Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000 * (cfg.BR/cfg.BR0) * R;
    A_g=2*pi*(r_m+t_g/2)*h_y1; Phi_g=Bg*A_g; By1=Phi_g/(pi*r_m^2);
    dcs=r_p*d_c;
    n_half=(h_c/dcs)*((2/sqrt(3))*(max(w_c,0)/dcs));
    L_m=2*pi*(r_m+c_c+t_b+w_c/2);
    n_eff=n_sym*(h_y1+0.08*h_m)/h_c*n_half; n_eff=min(n_eff, n_sym*n_half);
    Kf=n_eff*Bg*L_m;
    m_coil=8960*(n_sym*n_half)*(pi/4)*d_c^2*L_m;

    Ia=b*t^3/12; k_x=24*mat.E*Ia/L^3;
    m_ms =cfg.HOLLOW*mat.rho*a_s^2*b; m_int=cfg.HOLLOW*2*mat.rho*d_b*h_i*b; m_bld=8*mat.rho*L*t*b;
    m_eff=m_ms+m_coil+0.25*m_int+(1/3)*m_bld;
    f1=(1/(2*pi))*sqrt(k_x/m_eff);
    x_max=Kf*I_max/k_x;
    sigma=3*mat.E*t*(x_max/2)/L^2;
    VCMd=2*(r_m+t_g+t_y3); footY=a_s+2*(2*L+h_i);

    r.f1=f1; r.k_x=k_x; r.m_eff=m_eff; r.x_max=x_max; r.sigma=sigma;
    r.Kf=Kf; r.Bg=Bg; r.By1=By1; r.VCMd=VCMd; r.footY=footY; r.m_coil=m_coil;
    r.w_c=w_c; r.dcs=dcs; r.n=n_sym*n_half; r.R=R;
end

function [c,ceq]=nlcon(xv,mat,cfg)
    r=perf(xv,mat,cfg); g=xv(8); x_req=1e-3; f_req=100;
    c(1)=g*x_req - r.x_max;
    c(2)=g*f_req - r.f1;
    c(3)=cfg.KT*r.sigma - 0.3*mat.Sy;
    c(4)=r.By1 - cfg.BY1MAX;
    c(5)=r.VCMd - 20e-3;
    c(6)=r.footY - 200e-3;
    c(7)=r.dcs - r.w_c;
    c(8)=r.Bg - 0.8;
    ceq=[];
end
