%% material_study.m  - 재료/자기회로/경량화 민감도 (surrogate B_g 기준, d_c=0.5 고정)
%  "재료를 바꾸면 KPI를 충족할 수 있는가?"  baseline은 본문 최적화와 동일한 FEMM surrogate.
%  레버: (i) 유연기구 재료, (ii) 자석/요크 업그레이드 bg_scale(N52+FeCo~1.3x), (iii) stage 중공 fstage.
clear; clc; rng(7);
Cs=[199.0491,92.0382,-30.2822,-82.8183,-0.3287,1.0510,8.0321,-3.1713,-7.3583,4.5247]; % FEMM B_g surrogate [mT]
MAT.Al={71e9,2770,276e6}; MAT.Mg={45e9,1740,200e6}; MAT.Be={287e9,1850,240e6}; MAT.Ti={114e9,4430,880e6};

cases={ % 라벨                       재료  bg_scale fstage
 'Al6061 (기준)',                 'Al', 1.0, 1.0;
 'Mg(AZ31)로 교체',               'Mg', 1.0, 1.0;
 'Be로 교체',                     'Be', 1.0, 1.0;
 'Ti-6Al-4V로 교체',              'Ti', 1.0, 1.0;
 'Al + 자석/요크 (N52/FeCo)',     'Al', 1.3, 1.0;
 'Al + stage 중공화',             'Al', 1.0, 0.4;
 'Mg + 자석 + 중공',              'Mg', 1.3, 0.4;
 'Be + 자석 + 중공',              'Be', 1.3, 0.4};

% [t L b r_m h_m h_y1 t_g d_c gamma], d_c=0.5 고정
lb=[0.3,20,8,2.0,25,4,1.0,0.5,0.2]; ub=[1.5,45,18,6.5,25,16,3.0,0.5,3.0]; x0=[0.8,40,15,4,25,9,2.0,0.5,0.5];
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',200,...
    'MaxFunctionEvaluations',2000,'ConstraintTolerance',1e-9);

fprintf('=== 재료/자기회로/경량화 민감도 (surrogate, d_c=0.5mm) ===\n');
fprintf('%-26s%8s%9s%10s%9s%8s\n','케이스','gamma','f1[Hz]','stroke','m_eff','판정');
fprintf('%s\n',repmat('-',1,72));
for i=1:size(cases,1)
    mp=MAT.(cases{i,2}); mat.E=mp{1}; mat.rho=mp{2}; mat.Sy=mp{3};
    mat.bg=cases{i,3}; mat.fs=cases{i,4}; mat.Cs=Cs;
    best=0.2; bx=x0; starts=[x0; lb+rand(15,9).*(ub-lb)];
    for s=1:size(starts,1)
        try [xo,fv,ef]=fmincon(@(x)-x(9),starts(s,:),[],[],[],[],lb,ub,@(x)nlc(x,mat),opts);
            if ef>0,[c,~]=nlc(xo,mat); if all(c<=1e-6)&&-fv>best,best=-fv;bx=xo;end; end
        catch; end
    end
    r=perf(bx,mat); v='미충족'; if best>=1.0,v='충족 OK'; end
    fprintf('%-26s%8.3f%9.1f%9.3fmm%8.1fg%9s\n',cases{i,1},best,r.f1,r.x_max*1e3,r.m_eff*1e3,v);
end

function r=perf(xv,mat)
    t=xv(1)*1e-3;L=xv(2)*1e-3;b=xv(3)*1e-3;r_m=xv(4)*1e-3;h_m=xv(5)*1e-3;h_y1=xv(6)*1e-3;t_g=xv(7)*1e-3;d_c=xv(8)*1e-3;
    c_c=0.3e-3;t_b=0.3e-3;r_p=1.15;n_sym=2;I_max=2; a_s=15e-3;d_b=40e-3;h_i=5e-3; w_c=t_g-2*c_c;h_c=h_y1;
    rm=r_m*1e3;hy=h_y1*1e3;tg=t_g*1e3;
    Bg=mat.bg*(mat.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
    A_g=2*pi*(r_m+t_g/2)*h_y1; By1=Bg*A_g/(pi*r_m^2);
    dcs=r_p*d_c; n_half=(h_c/dcs)*((2/sqrt(3))*(max(w_c,0)/dcs)); L_m=2*pi*(r_m+c_c+t_b+w_c/2);
    n_eff=min(n_sym*(h_y1+0.08*h_m)/h_c*n_half,n_sym*n_half); Kf=n_eff*Bg*L_m;
    m_coil=8960*(n_sym*n_half)*(pi/4)*d_c^2*L_m;
    Ia=b*t^3/12; k_x=24*mat.E*Ia/L^3;
    m_eff=mat.fs*mat.rho*a_s^2*b+m_coil+mat.fs*0.25*2*mat.rho*d_b*h_i*b+(1/3)*8*mat.rho*L*t*b;
    f1=(1/(2*pi))*sqrt(k_x/m_eff); x_max=Kf*I_max/k_x; sigma=3*mat.E*t*(x_max/2)/L^2;
    r.f1=f1;r.x_max=x_max;r.sigma=sigma;r.By1=By1;r.Bg=Bg;r.w_c=w_c;r.dcs=dcs;r.m_eff=m_eff;
    r.VCMd=2*(r_m+t_g+1e-3); r.footY=a_s+2*(2*L+h_i);
end
function [c,ceq]=nlc(xv,mat)
    r=perf(xv,mat); g=xv(9);
    c=[g*1e-3-r.x_max; g*100-r.f1; r.sigma-0.3*mat.Sy; r.By1-1.8; r.VCMd-20e-3; r.footY-200e-3; r.dcs-r.w_c]; ceq=[];
end
