%% codesign_dc_material.m - d_c{0.3,0.5} x 실재료 8종 공동최적화 (leaf)
%  스펙(d_c>0.5) 준수 가능성을 확장 재료로 재판정
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2;
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',250,...
     'MaxFunctionEvaluations',2500,'ConstraintTolerance',1e-9);
DB={'Mg AZ31',45,200,1740; 'Al6061-T6',69,276,2700; 'Al7075-T6',72,503,2810;
    'Ti-6Al-4V',114,880,4430; 'BeCu C17200',128,1100,8250; '17-4PH H900',197,1170,7800;
    'Spring steel',200,1200,7850; 'Maraging300',190,2000,8100};

fprintf('=== d_c 0.3 vs 0.5 x 실재료 8종 (leaf, 피로0.2Sy, N52) ===\n');
fprintf('%-14s%14s%14s\n','재료','gamma(dc=0.3)','gamma(dc=0.5)');
fprintf('%s\n',repmat('-',1,44));
best5=-Inf; best5n='';
for k=1:size(DB,1)
  m.E=DB{k,2}*1e9; m.Sy=DB{k,3}*1e6; m.rho=DB{k,4};
  g3=solveleaf(m,0.3e-3,P,opts,14);
  g5=solveleaf(m,0.5e-3,P,opts,14);
  v3=tagv(g3); v5=tagv(g5);
  fprintf('%-14s%10.3f%s%10.3f%s\n', DB{k,1}, g3,v3, g5,v5);
  if g5>best5, best5=g5; best5n=DB{k,1}; end
end
fprintf('%s\n',repmat('-',1,44));
fprintf('d_c=0.5(스펙) 최고: %s, gamma=%.3f -> %s\n', best5n, best5, ...
   ternary(best5>=1.0,'스펙 충족 가능!','여전히 미충족'));

%% ===== 함수 =====
function s=tagv(g), if g>=1.0, s=' OK'; else, s='   '; end, end
function s=ternary(c,a,b), if c, s=a; else, s=b; end, end

function g=solveleaf(m,dc,P,opts,ns)
  lb=[0.3,15,5,15,5,2,1.2,0.1]; ub=[2.0,60,20,45,8,7,2.5,3]; x0=[0.6,35,10,25,7,4,1.5,0.7];
  bg=-Inf; S=[x0; lb+rand(ns,8).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(8),S(s,:),[],[],[],[],lb,ub,@(x)conL(x,m,dc,P),opts);
      if ef>0, c=conL(xo,m,dc,P); if all(c<=1e-6)&&(-fv)>bg, bg=-fv; end, end
    catch
    end
  end
  g=bg;
end
function v=vcm(rm,hy,tg,dc,P)
  c_c=0.3e-3;t_b=0.3e-3;rp=1.15;hm=40e-3;ty3=0.8e-3;nsym=2;
  Bg=(P.Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
  rmM=rm/1e3;hyM=hy/1e3;tgM=tg/1e3; wc=tgM-2*c_c; dcs=rp*dc; hc=hyM;
  n_half=(hc/dcs)*((2/sqrt(3))*(max(wc,0)/dcs)); L_m=2*pi*(rmM+c_c+t_b+wc/2);
  n_eff=min(nsym*(hyM+0.08*hm)/hc*n_half, nsym*n_half);
  v.Kf=n_eff*Bg*L_m; v.mcoil=8960*(nsym*n_half)*(pi/4)*dc^2*L_m;
  A_g=2*pi*(rmM+tgM/2)*hyM; v.By1=Bg*A_g/(pi*rmM^2); v.VCMd=2*(rmM+tgM+ty3); v.wc=wc; v.dcs=dcs;
end
function r=perfL(tm,Lm,bm,am,rm,hy,tg,E,rho,dc,P)
  t=tm/1e3;L=Lm/1e3;b=bm/1e3;a=am/1e3; v=vcm(rm,hy,tg,dc,P); db=40e-3;hi=5e-3;
  Keff=2*E*b*t^3/L^3; Mflex=rho*a^2*b+0.25*(2*rho*db*hi*b)+(1/3)*(8*rho*L*t*b);
  Meff=Mflex+v.mcoil; xs=v.Kf*P.I/Keff; f1=(1/(2*pi))*sqrt(Keff/Meff);
  r.x=xs;r.f1=f1;r.By1=v.By1;r.VCMd=v.VCMd;r.dcs=v.dcs;r.wc=v.wc;r.t=t;r.L=L;r.foot=a+2*(2*L+hi);
end
function [c,ceq]=conL(x,m,dc,P)
  r=perfL(x(1),x(2),x(3),x(4),x(5),x(6),x(7),m.E,m.rho,dc,P); g=x(8); xop=g*1e-3;
  peak=P.KT*3*m.E*r.t*(xop/2)/r.L^2;
  c=[g*1e-3-r.x; g*100-r.f1; peak-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2];
  ceq=[];
end
