%% opt_femm_n52_system.m
%  FEMM 실측 N52 surrogate(R^2=0.996)로 시스템 재최적화
%  d_c=0.5 고정, 솔리드(중공 없음), r_m<=8(Ø20 한계까지), K_t=1.7 피로 포함
%  VCM: h_m=40, t_y3=0.8 고정(DOE와 동일). 변수 x=[t,L,b,r_m,h_y1,t_g, gamma]
clear; clc; rng(7);
% N52 FEMM-적합 surrogate (mT), 입력 mm: [1,rm,hy,tg,rm^2,hy^2,tg^2,rm*hy,rm*tg,hy*tg]
Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
cfg.Cs=Cs; cfg.KT=1.7; cfg.DC=0.5e-3; cfg.AS=15e-3;   % a_s=15mm (경량화 트릭 없음)

opts=optimoptions('fmincon','Display','off','Algorithm','sqp',...
     'MaxIterations',300,'MaxFunctionEvaluations',3000,'ConstraintTolerance',1e-9);
% [ t    L    b    r_m  h_y1 t_g  gamma]
lb=[0.3, 20,  8,   5.0, 2.0, 1.2, 0.2];
ub=[1.5, 45,  18,  8.0, 7.0, 2.5, 3.0];
x0=[0.4, 28,  8,   8.0, 3.0, 1.2, 0.9];

mats={'Al6061',71e9,2770,276e6; 'Mg AZ31',45e9,1740,200e6};
fprintf('=== FEMM-N52 surrogate 시스템 재최적화 (d_c=0.5, 솔리드, K_t=1.7) ===\n\n');
fprintf('%-10s%8s%11s%9s%9s%8s%8s%8s%9s\n','재료','gamma','stroke[mm]','f1[Hz]','Bg[mT]','Kf','r_m','t_g','m_eff[g]');
fprintf('%s\n',repmat('-',1,80));
for k=1:2
  mat.name=mats{k,1}; mat.E=mats{k,2}; mat.rho=mats{k,3}; mat.Sy=mats{k,4};
  [bx,br]=solve_opt(lb,ub,x0,mat,cfg,opts,18);
  v=''; if bx(7)>=1.0, v=' <= 충족!'; end
  fprintf('%-10s%8.3f%11.3f%9.1f%9.0f%8.2f%8.2f%8.2f%9.1f%s\n', ...
     mat.name,bx(7),br.x_max*1e3,br.f1,br.Bg*1e3,br.Kf,bx(4),bx(6),br.m_eff*1e3,v);
  if k==2 || bx(7)>=1.0
    fprintf('   상세: t=%.2f L=%.1f b=%.1f mm, h_y1=%.2f | By1=%.2f T, n=%.0f turns, m_coil=%.1fg, Ø=%.1fmm\n',...
       bx(1),bx(2),bx(3),bx(5),br.By1,round(br.n),br.m_coil*1e3,br.VCMd*1e3);
    fprintf('   피로: peak(Kt)=%.1f MPa, SF=%.2f\n', cfg.KT*br.sigma*1e-6, 0.3*mat.Sy/(cfg.KT*br.sigma));
  end
end

%% ===== 함수 =====
function [bestx,bestr]=solve_opt(lb,ub,x0,mat,cfg,opts,nstart)
  bestg=-Inf; bestx=x0;
  starts=[x0; lb+rand(nstart,7).*(ub-lb)];
  for s=1:size(starts,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(7), starts(s,:),[],[],[],[],lb,ub,@(x)nlc(x,mat,cfg),opts);
      if ef>0, c=nlc(xo,mat,cfg); if all(c<=1e-6)&&(-fv)>bestg, bestg=-fv; bestx=xo; end, end
    catch
    end
  end
  bestr=perf(bestx,mat,cfg);
end

function r=perf(xv,mat,cfg)
  t=xv(1)*1e-3; L=xv(2)*1e-3; b=xv(3)*1e-3; r_m=xv(4)*1e-3; h_y1=xv(5)*1e-3; t_g=xv(6)*1e-3;
  d_c=cfg.DC; h_m=40e-3; t_y3=0.8e-3; c_c=0.3e-3; t_b=0.3e-3; r_p=1.15;
  n_sym=2; I_max=2; a_s=cfg.AS; d_b=40e-3; h_i=5e-3; w_c=t_g-2*c_c; h_c=h_y1;
  rm=r_m*1e3; hy=h_y1*1e3; tg=t_g*1e3;
  Bg=(cfg.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;   % N52 FEMM surrogate
  A_g=2*pi*(r_m+t_g/2)*h_y1; Phi_g=Bg*A_g; By1=Phi_g/(pi*r_m^2);
  dcs=r_p*d_c; n_half=(h_c/dcs)*((2/sqrt(3))*(max(w_c,0)/dcs));
  L_m=2*pi*(r_m+c_c+t_b+w_c/2);
  n_eff=min(n_sym*(h_y1+0.08*h_m)/h_c*n_half, n_sym*n_half);
  Kf=n_eff*Bg*L_m; m_coil=8960*(n_sym*n_half)*(pi/4)*d_c^2*L_m;
  Ia=b*t^3/12; k_x=24*mat.E*Ia/L^3;
  m_ms=mat.rho*a_s^2*b; m_int=2*mat.rho*d_b*h_i*b; m_bld=8*mat.rho*L*t*b;  % solid
  m_eff=m_ms+m_coil+0.25*m_int+(1/3)*m_bld;
  f1=(1/(2*pi))*sqrt(k_x/m_eff); x_max=Kf*I_max/k_x;
  sigma=3*mat.E*t*(x_max/2)/L^2; VCMd=2*(r_m+t_g+t_y3); footY=a_s+2*(2*L+h_i);
  r.f1=f1;r.k_x=k_x;r.m_eff=m_eff;r.x_max=x_max;r.sigma=sigma;r.Kf=Kf;r.Bg=Bg;
  r.By1=By1;r.VCMd=VCMd;r.footY=footY;r.m_coil=m_coil;r.w_c=w_c;r.dcs=dcs;r.n=n_sym*n_half;
end

function [c,ceq]=nlc(xv,mat,cfg)
  r=perf(xv,mat,cfg); g=xv(7);
  c(1)=g*1e-3 - r.x_max;
  c(2)=g*100  - r.f1;
  c(3)=cfg.KT*r.sigma - 0.3*mat.Sy;
  c(4)=r.By1 - 2.2;
  c(5)=r.VCMd - 20e-3;
  c(6)=r.footY - 200e-3;
  c(7)=r.dcs - r.w_c;
  c(8)=r.Bg - 1.1;
  ceq=[];
end
