%% opt_design.m
%  VCM 구동 Double Compound Parallelogram 스테이지 — 시스템 최적화
%  목적: 스펙 충족률 gamma 최대화 (epigraph): stroke>=gamma*1mm, f1>=gamma*100Hz
%        gamma>=1 이면 두 스펙 동시 충족, <1 이면 최선의 충족 비율(미충족 정도).
%  설계변수 x=[t,L,b,r_m,h_m,h_y1,t_g,d_c]  (mm),  + gamma
%  나머지 제약: 응력<=0.3Sy, 요크 B<=1.8T, VCM직경<=20mm, footprintY<=200mm, 코일이 공극에 들어감
clear; clc; rng(7);

mats = { 'Al6061-T6',    71e9,   2770, 276e6;
         'Spring-steel', 200e9,  7850, 1200e6;
         'BeCu C17200',  128e9,  8250, 1100e6 };

% [ t    L    b    r_m  h_m  h_y1 t_g  d_c  gamma]   (h_m 고정25, r_m/h_y1/t_g는 surrogate 격자 박스 내)
lb = [0.3, 20,  8,   2.0, 25,  4,   1.0, 0.5, 0.2];
ub = [1.5, 45,  18,  6.5, 25,  16,  3.0, 1.2, 3.0];
x0 = [0.8, 40,  15,  4,   25,  9,   2.0, 0.6, 0.5];

opts = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
        'MaxIterations',250,'MaxFunctionEvaluations',2500, ...
        'StepTolerance',1e-11,'ConstraintTolerance',1e-9);

%% ===== 메인: d_c>=0.5mm (스펙대로) =====
fprintf('=== [A] 스펙대로 d_c>=0.5mm : 스펙 충족률 gamma 최대화 ===\n');
fprintf('%-13s%8s%10s%9s%9s%8s%8s%8s%9s\n', ...
    'Material','gamma','stroke[mm]','f1[Hz]','sig[MPa]','SF','Kf[N/A]','d_c[mm]','m_eff[g]');
fprintf('%s\n',repmat('-',1,84));
BX=cell(3,1); BR=cell(3,1);
for k=1:3
    mat.name=mats{k,1}; mat.E=mats{k,2}; mat.rho=mats{k,3}; mat.Sy=mats{k,4};
    [bx,br]=solve_opt(lb,ub,x0,mat,opts,10); BX{k}=bx; BR{k}=br;
    fprintf('%-13s%8.3f%10.3f%9.1f%9.1f%8.2f%8.3f%8.2f%9.1f\n', ...
        mat.name, bx(9), br.x_max*1e3, br.f1, br.sigma*1e-6, br.SF, br.Kf, bx(8), br.m_eff*1e3);
end
fprintf('\n--- 최적해 상세 (물리 정합성 확인) ---\n');
fprintf('%-13s%7s%7s%7s%7s%7s%7s%8s%8s%8s%9s%9s\n', ...
    'Material','t','L','b','r_m','h_m','h_y1','t_g','Bg[mT]','By1[T]','n(turns)','mcoil[g]');
for k=1:3
    bx=BX{k}; br=BR{k};
    fprintf('%-13s%7.2f%7.1f%7.1f%7.2f%7.1f%7.1f%8.2f%8.1f%8.2f%9.0f%9.1f\n', ...
        mats{k,1}, bx(1),bx(2),bx(3),bx(4),bx(5),bx(6),bx(7), br.Bg*1e3, br.By1, round(br.n), br.m_coil*1e3);
end

%% ===== d_c 민감도 (Al6061) : 어느 선경에서 gamma>=1 회복? =====
fprintf('\n=== [B] d_c 민감도 (Al6061, d_c 고정) : 스펙 충족 임계 선경 ===\n');
fprintf('%-9s%8s%10s%9s%9s\n','d_c[mm]','gamma','stroke[mm]','f1[Hz]','판정');
fprintf('%s\n',repmat('-',1,46));
mat.name='Al6061-T6'; mat.E=71e9; mat.rho=2770; mat.Sy=276e6;
for dc=[0.50 0.35 0.30 0.26 0.22]
    lb2=lb; ub2=ub; lb2(8)=dc; ub2(8)=dc; x0b=x0; x0b(8)=dc;
    [bx,br]=solve_opt(lb2,ub2,x0b,mat,opts,8);
    verdict='미충족'; if bx(9)>=1.0, verdict='충족 OK'; end
    fprintf('%-9.2f%8.3f%10.3f%9.1f%9s\n', dc, bx(9), br.x_max*1e3, br.f1, verdict);
end

fprintf('\n주: gamma>=1 => stroke>=1mm & f1>=100Hz 동시 충족. gamma<1 => 그 비율만 달성(미충족).\n');

%% ================== 함수 ==================
function [bestx,bestr]=solve_opt(lb,ub,x0,mat,opts,nstart)
    bestg=-Inf; bestx=x0;
    starts=[x0; lb+rand(nstart,9).*(ub-lb)];
    for s=1:size(starts,1)
        try
            [xo,fv,ef]=fmincon(@(x)-x(9), starts(s,:), [],[],[],[], lb,ub, @(x)nlcon(x,mat), opts);
            if ef>0
                [c,~]=nlcon(xo,mat);
                if all(c<=1e-6) && (-fv)>bestg, bestg=-fv; bestx=xo; end
            end
        catch
        end
    end
    bestr=perf(bestx,mat);
end

function r=perf(xv,mat)
    mu0=4*pi*1e-7;
    t=xv(1)*1e-3; L=xv(2)*1e-3; b=xv(3)*1e-3;
    r_m=xv(4)*1e-3; h_m=xv(5)*1e-3; h_y1=xv(6)*1e-3; t_g=xv(7)*1e-3; d_c=xv(8)*1e-3;
    t_y3=1e-3; h_y2=3e-3; c_c=0.3e-3; t_b=0.3e-3; r_p=1.15;
    Br=1.17; Hc=867e3; mu_s=1400; n_sym=2; I_max=2;
    a_s=15e-3; d_b=40e-3; h_i=5e-3; x_req=1e-3;
    w_c=t_g-2*c_c; h_c=h_y1;

    mu_r=Br/(mu0*Hc);
    P_m =pi*mu0*mu_r*r_m^2/h_m;
    P_y1=pi*mu0*mu_s*r_m^2/h_y1;
    P_y2=2*pi*mu0*mu_s*h_y2/(log(r_m+t_g+t_y3/2)-log(r_m/2));
    P_y3=pi*mu0*mu_s*((r_m+t_g+t_y3)^2-(r_m+t_g)^2)/(h_y1+h_m);
    P_g =2*pi*mu0*h_y1/(log(r_m+t_g)-log(r_m));
    P_l =0.52*(2*pi*mu0*r_m);
    % B_g: FEMM-적합 surrogate (bg_surrogate.py, R^2=0.99). 입력 mm -> 출력 mT -> T
    Cs=[199.0491,92.0382,-30.2822,-82.8183,-0.3287,1.0510,8.0321,-3.1713,-7.3583,4.5247];
    rm=r_m*1e3; hy=h_y1*1e3; tg=t_g*1e3;
    Bg=(Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
    A_g=2*pi*(r_m+t_g/2)*h_y1; Phi_g=Bg*A_g; By1=Phi_g/(pi*r_m^2);
    dcs=r_p*d_c;
    % 최적화 평활화 위해 권선수 연속(최종 설계에서만 반올림). floor 제거.
    n_half=(h_c/dcs)*((2/sqrt(3))*(max(w_c,0)/dcs));
    L_m=2*pi*(r_m+c_c+t_b+w_c/2);
    n_eff=n_sym*(h_y1+0.08*h_m)/h_c*n_half;
    n_eff=min(n_eff, n_sym*n_half);        % 유효권선수 <= 실제권선수 (물리 가드)
    Kf=n_eff*Bg*L_m;
    m_coil=8960*(n_sym*n_half)*(pi/4)*d_c^2*L_m;

    Ia=b*t^3/12; k_x=24*mat.E*Ia/L^3;
    m_ms=mat.rho*a_s^2*b; m_int=2*mat.rho*d_b*h_i*b; m_bld=8*mat.rho*L*t*b;
    m_eff=m_ms+m_coil+0.25*m_int+(1/3)*m_bld;
    f1=(1/(2*pi))*sqrt(k_x/m_eff);
    F_max=Kf*I_max; x_max=F_max/k_x;
    sigma=3*mat.E*t*(x_max/2)/L^2; SF=mat.Sy/sigma;
    VCMd=2*(r_m+t_g+t_y3); footY=a_s+2*(2*L+h_i);

    r.f1=f1; r.k_x=k_x; r.m_eff=m_eff; r.x_max=x_max; r.sigma=sigma; r.SF=SF;
    r.Kf=Kf; r.Bg=Bg; r.By1=By1; r.VCMd=VCMd; r.footY=footY; r.m_coil=m_coil;
    r.w_c=w_c; r.dcs=dcs; r.n=n_sym*n_half;
end

function [c,ceq]=nlcon(xv,mat)
    r=perf(xv,mat); g=xv(9); x_req=1e-3; f_req=100;
    c(1)=g*x_req - r.x_max;     % stroke >= gamma*1mm
    c(2)=g*f_req - r.f1;        % f1 >= gamma*100Hz
    c(3)=r.sigma - 0.3*mat.Sy;  % 응력
    c(4)=r.By1 - 1.8;           % 포화
    c(5)=r.VCMd - 20e-3;        % VCM 직경
    c(6)=r.footY - 200e-3;      % footprint
    c(7)=r.dcs - r.w_c;         % 코일 폭 w_c >= d_c* (소선 최소 1가닥)
    c(8)=r.Bg - 0.6;            % 물리 sanity 상한 (surrogate가 정확하므로 가드 거의 불필요)
    ceq=[];
end
