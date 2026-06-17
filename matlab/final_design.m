%% final_design.m
%  최종 설계안 — Al6061 기준, d_c 해석별 완성 설계 시트
%  d_c=0.5(소선해석/literal), 0.4(절충), 0.3(보어해석/완전충족)
clear; clc; rng(7);

mat.name='Al6061-T6'; mat.E=71e9; mat.rho=2770; mat.Sy=276e6;
% [ t    L    b    r_m  h_m  h_y1 t_g  d_c  gamma]
lb=[0.3,20, 8, 2.0,25,4, 1.0,0.5,0.2]; ub=[1.5,45,18,6.5,25,16,3.0,1.2,3.0]; % h_m25고정, surrogate박스
x0=[0.8,40,15,4,25,9,2.0,0.6,0.5];
opts=optimoptions('fmincon','Display','off','Algorithm','sqp', ...
     'MaxIterations',250,'MaxFunctionEvaluations',2500,'StepTolerance',1e-11,'ConstraintTolerance',1e-9);

for dc=[0.50 0.40 0.30]
    lb2=lb; ub2=ub; lb2(8)=dc; ub2(8)=dc; x0b=x0; x0b(8)=dc;
    bestg=-Inf; bx=x0b; starts=[x0b; lb2+rand(15,9).*(ub2-lb2)];
    for s=1:size(starts,1)
        try
            [xo,fv,ef]=fmincon(@(x)-x(9),starts(s,:),[],[],[],[],lb2,ub2,@(x)nlcon(x,mat),opts);
            if ef>0, [c,~]=nlcon(xo,mat); if all(c<=1e-6)&&(-fv)>bestg, bestg=-fv; bx=xo; end; end
        catch; end
    end
    r=perf(bx,mat); I1=r.k_x*1e-3/r.Kf; P1=I1^2*r.Rc; P2=2^2*r.Rc;
    fprintf('================ d_c = %.2f mm =================\n', dc);
    fprintf(' [유연기구]  t=%.2f mm  L=%.1f mm  b=%.1f mm\n', bx(1),bx(2),bx(3));
    fprintf(' [VCM]      r_m=%.2f  h_m=%.1f  h_y1=%.1f  t_g=%.2f mm  (직경 %.1f mm)\n', ...
             bx(4),bx(5),bx(6),bx(7), r.VCMd*1e3);
    fprintf(' 강성 k_x   = %.0f N/m\n', r.k_x);
    fprintf(' 공진 f1    = %.1f Hz   (요구 >100)\n', r.f1);
    fprintf(' 최대변위   = %.3f mm @2A  (요구 >1)\n', r.x_max*1e3);
    fprintf(' 1mm 전류   = %.2f A    (요구 <=2)\n', I1);
    fprintf(' blade 응력 = %.1f MPa  (허용 0.3Sy=%.0f), SF=%.2f\n', r.sigma*1e-6, 0.3*mat.Sy*1e-6, r.SF);
    fprintf(' 힘상수 Kf  = %.3f N/A   Bg=%.0f mT  By1=%.2f T  권선 n=%d\n', r.Kf, r.Bg*1e3, r.By1, round(r.n));
    fprintf(' 질량       : 코일 %.1f g, m_eff %.1f g\n', r.m_coil*1e3, r.m_eff*1e3);
    fprintf(' 발열       : 1mm유지 %.2f W, 2A최대 %.2f W\n', P1, P2);
    fprintf(' footprint  : Y=%.0f mm (요구 <=200), 높이 b=%.0f, VCM직경=%.1f (요구<=20)\n\n', ...
             r.footY*1e3, bx(3), r.VCMd*1e3);
end

%% ===== 모델 (opt_design.m과 동일) =====
function r=perf(xv,mat)
    mu0=4*pi*1e-7;
    t=xv(1)*1e-3;L=xv(2)*1e-3;b=xv(3)*1e-3;r_m=xv(4)*1e-3;h_m=xv(5)*1e-3;h_y1=xv(6)*1e-3;t_g=xv(7)*1e-3;d_c=xv(8)*1e-3;
    t_y3=1e-3;h_y2=3e-3;c_c=0.3e-3;t_b=0.3e-3;r_p=1.15;Br=1.17;Hc=867e3;mu_s=1400;n_sym=2;I_max=2;
    a_s=15e-3;d_b=40e-3;h_i=5e-3; w_c=t_g-2*c_c;h_c=h_y1;
    mu_r=Br/(mu0*Hc);
    P_m=pi*mu0*mu_r*r_m^2/h_m; P_y1=pi*mu0*mu_s*r_m^2/h_y1;
    P_y2=2*pi*mu0*mu_s*h_y2/(log(r_m+t_g+t_y3/2)-log(r_m/2));
    P_y3=pi*mu0*mu_s*((r_m+t_g+t_y3)^2-(r_m+t_g)^2)/(h_y1+h_m);
    P_g=2*pi*mu0*h_y1/(log(r_m+t_g)-log(r_m)); P_l=0.52*(2*pi*mu0*r_m);
    Cs=[199.0491,92.0382,-30.2822,-82.8183,-0.3287,1.0510,8.0321,-3.1713,-7.3583,4.5247]; % FEMM surrogate
    rm=r_m*1e3; hy=h_y1*1e3; tg=t_g*1e3;
    Bg=(Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
    A_g=2*pi*(r_m+t_g/2)*h_y1; Phi_g=Bg*A_g; By1=Phi_g/(pi*r_m^2);
    dcs=r_p*d_c; n_half=(h_c/dcs)*((2/sqrt(3))*(max(w_c,0)/dcs)); L_m=2*pi*(r_m+c_c+t_b+w_c/2);
    n_eff=n_sym*(h_y1+0.08*h_m)/h_c*n_half; n_eff=min(n_eff,n_sym*n_half); Kf=n_eff*Bg*L_m;
    m_coil=8960*(n_sym*n_half)*(pi/4)*d_c^2*L_m; Rc=1.68e-8*(n_sym*n_half)*L_m/((pi/4)*d_c^2);
    Ia=b*t^3/12; k_x=24*mat.E*Ia/L^3;
    m_ms=mat.rho*a_s^2*b; m_int=2*mat.rho*d_b*h_i*b; m_bld=8*mat.rho*L*t*b;
    m_eff=m_ms+m_coil+0.25*m_int+(1/3)*m_bld; f1=(1/(2*pi))*sqrt(k_x/m_eff);
    F_max=Kf*I_max; x_max=F_max/k_x; sigma=3*mat.E*t*(x_max/2)/L^2; SF=mat.Sy/sigma;
    VCMd=2*(r_m+t_g+t_y3); footY=a_s+2*(2*L+h_i);
    r.f1=f1;r.k_x=k_x;r.m_eff=m_eff;r.x_max=x_max;r.sigma=sigma;r.SF=SF;r.Kf=Kf;r.Bg=Bg;r.By1=By1;
    r.VCMd=VCMd;r.footY=footY;r.m_coil=m_coil;r.Rc=Rc;r.w_c=w_c;r.dcs=dcs;r.n=n_sym*n_half;
end
function [c,ceq]=nlcon(xv,mat)
    r=perf(xv,mat); g=xv(9);
    c(1)=g*1e-3-r.x_max; c(2)=g*100-r.f1; c(3)=r.sigma-0.3*mat.Sy; c(4)=r.By1-1.8;
    c(5)=r.VCMd-20e-3; c(6)=r.footY-200e-3; c(7)=r.dcs-r.w_c; c(8)=r.Bg-0.6; ceq=[]; % Bg는 FEMM surrogate(가드 거의 불필요)
end
