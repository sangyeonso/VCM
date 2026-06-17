%% codesign_mono.m  -  유연기구 + VCM 일괄(monolithic) 공동최적화
%  트랙: hinge{leaf,notch} x d_c{0.3,0.5}, 재료 내부루프, 목적 max gamma
%  시스템 커플링: x=Kf*I/Keff,  f1=(1/2pi)sqrt(Keff/(M_flex+m_coil))
%  제약: x>=g*1mm, f1>=g*100Hz, 응력(피로 0.2Sy, Kt), VCM Ø<=20, foot<=200, 코일끼움, 포화
%  VCM B_g = FEMM-적합 N52 surrogate (R^2=0.996)
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2;          % 피로 허용 0.2*Sy, leaf 응력집중 Kt=1.7
mats={'Al',72e9,270e6,2770;'Mg',45e9,200e6,1740;'BeCu',128e9,1100e6,8250;'Spring',200e9,1200e6,7850};
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',250,...
     'MaxFunctionEvaluations',2500,'ConstraintTolerance',1e-9);

fprintf('=== 일괄 공동최적화 (피로 0.2Sy, Kt=1.7, N52 FEMM surrogate) ===\n');
fprintf('%-7s%6s%8s%8s%11s%9s%8s%9s%9s\n','hinge','d_c','재료','gamma','stroke[mm]','f1[Hz]','Kf','Keff[N/mm]','mc[g]');
fprintf('%s\n',repmat('-',1,76));
for hinge={'leaf','notch'}
 for dc=[0.3e-3 0.5e-3]
  bestg=-Inf; bestrow='';
  for k=1:4
    m.E=mats{k,2}; m.Sy=mats{k,3}; m.rho=mats{k,4};
    if strcmp(hinge{1},'leaf')
      lb=[0.3,15,5,15,5,2,1.2,0.1]; ub=[2.0,60,20,45,8,7,2.5,3]; x0=[0.6,35,10,25,7,4,1.5,0.7];
      [bx,fv]=solveMS(@(x)perfL(x,m,dc,P),@(x)conL(x,m,dc,P),lb,ub,x0,opts,12);
      r=perfL(bx,m,dc,P); g=bx(8);
    else
      lb=[0.3,0.2,5,10,15,5,2,1.2,0.1]; ub=[5,1.5,20,75,45,8,7,2.5,3]; x0=[2,0.5,10,40,25,7,4,1.5,0.7];
      [bx,fv]=solveMS(@(x)perfN(x,m,dc,P),@(x)conN(x,m,dc,P),lb,ub,x0,opts,12);
      r=perfN(bx,m,dc,P); g=bx(9);
    end
    if g>bestg, bestg=g;
      bestrow=sprintf('%-7s%6.1f%8s%8.3f%11.3f%9.1f%8.2f%9.0f%9.1f', ...
        hinge{1}, dc*1e3, mats{k,1}, g, r.x*1e3, r.f1, r.Kf, r.Keff/1e3, r.mcoil*1e3);
    end
  end
  vd=''; if bestg>=1.0, vd=' <= 충족'; end
  fprintf('%s%s\n', bestrow, vd);
 end
end
fprintf('\n주: 각 행은 4재료 중 best. gamma>=1 이면 1mm&100Hz 동시충족(피로+Kt+커플링 포함).\n');

%% ===== solver =====
function [bx,bf]=solveMS(pf,cf,lb,ub,x0,opts,ns)
  n=numel(lb); bg=-Inf; bx=x0;
  S=[x0; lb+rand(ns,n).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(end),S(s,:),[],[],[],[],lb,ub,cf,opts);
      if ef>0, c=cf(xo); if all(c<=1e-6)&&(-fv)>bg, bg=-fv; bx=xo; end, end
    catch
    end
  end
  bf=bg;
end

%% ===== shared VCM (N52 FEMM surrogate) =====
function v=vcm(rm,hy,tg,dc,P)
  c_c=0.3e-3;t_b=0.3e-3;rp=1.15;hm=40e-3;ty3=0.8e-3;nsym=2;
  Bg=(P.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;   % rm,hy,tg [mm]->T
  rmM=rm/1e3;hyM=hy/1e3;tgM=tg/1e3;
  wc=tgM-2*c_c; dcs=rp*dc; hc=hyM;
  n_half=(hc/dcs)*((2/sqrt(3))*(max(wc,0)/dcs));
  L_m=2*pi*(rmM+c_c+t_b+wc/2);
  n_eff=min(nsym*(hyM+0.08*hm)/hc*n_half, nsym*n_half);
  Kf=n_eff*Bg*L_m; mcoil=8960*(nsym*n_half)*(pi/4)*dc^2*L_m;
  A_g=2*pi*(rmM+tgM/2)*hyM; By1=Bg*A_g/(pi*rmM^2);
  v.Bg=Bg;v.Kf=Kf;v.mcoil=mcoil;v.By1=By1;v.VCMd=2*(rmM+tgM+ty3);v.wc=wc;v.dcs=dcs;
end

%% ===== leaf flexure =====
function r=perfL(x,m,dc,P)
  t=x(1)/1e3;L=x(2)/1e3;b=x(3)/1e3;a=x(4)/1e3; v=vcm(x(5),x(6),x(7),dc,P);
  Keff=2*m.E*b*t^3/L^3;
  db=40e-3;hi=5e-3;
  M_flex=m.rho*a^2*b + 0.25*(2*m.rho*db*hi*b) + (1/3)*(8*m.rho*L*t*b);
  Meff=M_flex+v.mcoil; xs=v.Kf*P.I/Keff; f1=(1/(2*pi))*sqrt(Keff/Meff);
  r.Keff=Keff;r.Meff=Meff;r.x=xs;r.f1=f1;r.Kf=v.Kf;r.mcoil=v.mcoil;
  r.By1=v.By1;r.VCMd=v.VCMd;r.dcs=v.dcs;r.wc=v.wc;r.t=t;r.L=L;
  r.foot=a+2*(2*L+hi);
end
function [c,ceq]=conL(x,m,dc,P)
  r=perfL(x,m,dc,P); g=x(8); xop=g*1e-3;
  peak=P.KT*3*m.E*r.t*(xop/2)/r.L^2;
  c=[g*1e-3-r.x; g*100-r.f1; peak-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2];
  ceq=[];
end

%% ===== notch flexure (lecture) =====
function r=perfN(x,m,dc,P)
  R=x(1)/1e3;t=x(2)/1e3;b=x(3)/1e3;l=x(4)/1e3;a=x(5)/1e3; v=vcm(x(6),x(7),x(8),dc,P);
  Rmm=x(1);tmm=x(2);
  K=0.565*(tmm/Rmm)+0.166; Kt=(2.7*tmm+5.4*Rmm)/(8*Rmm+tmm)+0.325;
  k_th=2*m.E*b*t^2.5/(9*pi*sqrt(R)); Keff=4*k_th/l^2;
  thmax=(4*K/Kt)*(R/(m.E*t))*(P.SIG*m.Sy); xmax=2*l*thmax;
  M_flex=a^2*b*m.rho + 0.5*(a*3e-3*b*m.rho) + (2/3)*(l*3e-3*b*m.rho);
  Meff=M_flex+v.mcoil; xs=v.Kf*P.I/Keff; f1=(1/(2*pi))*sqrt(Keff/Meff);
  r.Keff=Keff;r.Meff=Meff;r.x=xs;r.f1=f1;r.xmax=xmax;r.Kf=v.Kf;r.mcoil=v.mcoil;
  r.By1=v.By1;r.VCMd=v.VCMd;r.dcs=v.dcs;r.wc=v.wc;r.R=x(1);r.t=x(2);
  r.foot=a+2*l;
end
function [c,ceq]=conN(x,m,dc,P)
  r=perfN(x,m,dc,P); g=x(9);
  c=[g*1e-3-r.x; g*100-r.f1; g*1e-3-r.xmax; r.t-r.R; r.R-5*r.t; ...
     r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2];
  ceq=[];
end
