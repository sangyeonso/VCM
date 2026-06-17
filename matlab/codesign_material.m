%% codesign_material.m - 재료를 최적화 파라미터로 (leaf, d_c=0.3)
%  (1) 확장 실재료 8종 각각 공동최적화 -> gamma 랭킹 + 물성지수
%  (2) 물성 (E,Sy,rho) 연속 최적화 -> 이상적 물성 타깃 (unobtainium 주의)
clear; clc; rng(7);
P.Cs=[489.4659,164.0726,-123.9167,-262.7587,-0.7897,7.3889,35.6957,-8.9347,-21.3500,28.9057];
P.I=2; P.KT=1.7; P.SIG=0.2; dc=0.3e-3;
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',250,...
     'MaxFunctionEvaluations',2500,'ConstraintTolerance',1e-9);

%% (1) 확장 실재료 DB  [name, E(GPa), Sy(MPa), rho]
DB={'Mg AZ31',45,200,1740; 'Al6061-T6',69,276,2700; 'Al7075-T6',72,503,2810;
    'Ti-6Al-4V',114,880,4430; 'BeCu C17200',128,1100,8250; '17-4PH H900',197,1170,7800;
    'Spring steel',200,1200,7850; 'Maraging300',190,2000,8100};
fprintf('=== (1) 실재료 8종 공동최적화 랭킹 (leaf, d_c=0.3) ===\n');
fprintf('%-14s%7s%7s%7s%9s%9s%9s%9s\n','재료','E','Sy','rho','Sy/E','sqrtE/rho','gamma','f1[Hz]');
fprintf('%s\n',repmat('-',1,74));
res=[];
for k=1:size(DB,1)
  m.E=DB{k,2}*1e9; m.Sy=DB{k,3}*1e6; m.rho=DB{k,4};
  [bx,r]=solveleaf(m,dc,P,opts,14);
  I1=m.Sy/m.E; I2=sqrt(m.E)/m.rho;            % 물성지수: 행정(Sy/E), 비강성(sqrtE/rho)
  fprintf('%-14s%7.0f%7.0f%7.0f%9.4f%9.4f%9.3f%9.1f\n',...
     DB{k,1},DB{k,2},DB{k,3},DB{k,4}, I1, I2*1e3, bx(end), r.f1);
  res(k)=bx(end);
end
[~,ix]=max(res); fprintf('-> 최고: %s (gamma=%.3f)\n', DB{ix,1}, res(ix));

%% (2) 물성 연속 최적화 : (E,Sy,rho)도 변수
fprintf('\n=== (2) 물성 연속 최적화 (E,Sy,rho 변수화) -> 이상적 타깃 ===\n');
% vars [t,L,b,a, rm,hy,tg, E(GPa),Sy(MPa),rho, g]
lb=[0.3,15,5,15, 5,2,1.2,  40,150,1700, 0.1];
ub=[2.0,60,20,45, 8,7,2.5, 210,2000,8300, 3];
x0=[0.6,35,10,25, 7,4,1.5,  70,300,2700, 0.7];
bg=-Inf; bx=x0; S=[x0; lb+rand(20,11).*(ub-lb)];
for s=1:size(S,1)
  try
    [xo,fv,ef]=fmincon(@(x)-x(11),S(s,:),[],[],[],[],lb,ub,@(x)conC(x,dc,P),opts);
    if ef>0, c=conC(xo,dc,P); if all(c<=1e-6)&&(-fv)>bg, bg=-fv; bx=xo; end, end
  catch
  end
end
fprintf('이상적 물성 타깃: E=%.0f GPa, Sy=%.0f MPa, rho=%.0f  -> gamma=%.3f\n',bx(8),bx(9),bx(10),bx(11));
fprintf('  지수: Sy/E=%.4f, sqrtE/rho=%.4f\n', (bx(9)*1e6)/(bx(8)*1e9), sqrt(bx(8)*1e9)/bx(10)*1e3);
fprintf('  해석: rho는 %s, Sy는 %s 으로 향함 (가볍고 강한 쪽). 실재료엔 없는 조합 가능성.\n',...
  ternary(bx(10)<2500,'최소(경량)','중간'), ternary(bx(9)>1500,'최대(고강도)','중간'));

%% ===== 함수 =====
function s=ternary(c,a,b), if c, s=a; else, s=b; end, end

function [bx,r]=solveleaf(m,dc,P,opts,ns)
  lb=[0.3,15,5,15,5,2,1.2,0.1]; ub=[2.0,60,20,45,8,7,2.5,3]; x0=[0.6,35,10,25,7,4,1.5,0.7];
  bg=-Inf; bx=x0; S=[x0; lb+rand(ns,8).*(ub-lb)];
  for s=1:size(S,1)
    try
      [xo,fv,ef]=fmincon(@(x)-x(8),S(s,:),[],[],[],[],lb,ub,@(x)conL(x,m,dc,P),opts);
      if ef>0, c=conL(xo,m,dc,P); if all(c<=1e-6)&&(-fv)>bg, bg=-fv; bx=xo; end, end
    catch
    end
  end
  r=perfL(bx(1),bx(2),bx(3),bx(4),bx(5),bx(6),bx(7),m.E,m.rho,dc,P);
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
  Keff=2*E*b*t^3/L^3;
  Mflex=rho*a^2*b + 0.25*(2*rho*db*hi*b) + (1/3)*(8*rho*L*t*b);
  Meff=Mflex+v.mcoil; xs=v.Kf*P.I/Keff; f1=(1/(2*pi))*sqrt(Keff/Meff);
  r.Keff=Keff;r.Meff=Meff;r.x=xs;r.f1=f1;r.Kf=v.Kf;r.mcoil=v.mcoil;
  r.By1=v.By1;r.VCMd=v.VCMd;r.dcs=v.dcs;r.wc=v.wc;r.t=t;r.L=L;r.foot=a+2*(2*L+hi);
end
function [c,ceq]=conL(x,m,dc,P)
  r=perfL(x(1),x(2),x(3),x(4),x(5),x(6),x(7),m.E,m.rho,dc,P); g=x(8); xop=g*1e-3;
  peak=P.KT*3*m.E*r.t*(xop/2)/r.L^2;
  c=[g*1e-3-r.x; g*100-r.f1; peak-P.SIG*m.Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2];
  ceq=[];
end
function [c,ceq]=conC(x,dc,P)   % 연속 물성: E=x(8)GPa,Sy=x(9)MPa,rho=x(10)
  E=x(8)*1e9; Sy=x(9)*1e6; rho=x(10);
  r=perfL(x(1),x(2),x(3),x(4),x(5),x(6),x(7),E,rho,dc,P); g=x(11); xop=g*1e-3;
  peak=P.KT*3*E*r.t*(xop/2)/r.L^2;
  c=[g*1e-3-r.x; g*100-r.f1; peak-P.SIG*Sy; r.VCMd-20e-3; r.foot-200e-3; r.dcs-r.wc; r.By1-2.2];
  ceq=[];
end
