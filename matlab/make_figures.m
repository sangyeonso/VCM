%% make_figures.m  - 보고서용 그림 생성 (figure 폴더에 PNG 저장)
%  라벨은 영어(폰트 안정), 캡션은 보고서에서 한글로.
clear; clc; close all;
FIGDIR = fullfile('..','figure');
Cs=[199.0491,92.0382,-30.2822,-82.8183,-0.3287,1.0510,8.0321,-3.1713,-7.3583,4.5247]; % B_g surrogate

%% ===== Fig 1. 재료별 스펙 충족률 gamma (d_c>=0.5) =====
mats={'Al6061',71e9,2770,276e6; 'Spring steel',200e9,7850,1200e6; 'BeCu',128e9,8250,1100e6};
lb=[0.3,20,8,2.0,25,4,1.0,0.5,0.2]; ub=[1.5,45,18,6.5,25,16,3.0,1.2,3.0];
x0=[0.8,40,15,4,25,9,2.0,0.6,0.5];
opts=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',200,...
    'MaxFunctionEvaluations',2000,'ConstraintTolerance',1e-9);
rng(7); gam=zeros(3,1);
for k=1:3
    mat.E=mats{k,2}; mat.rho=mats{k,3}; mat.Sy=mats{k,4};
    best=0; starts=[x0; lb+rand(25,9).*(ub-lb)];
    for s=1:size(starts,1)
        try [xo,fv,ef]=fmincon(@(x)-x(9),starts(s,:),[],[],[],[],lb,ub,@(x)nlc(x,mat,Cs),opts);
            if ef>0,[c,~]=nlc(xo,mat,Cs); if all(c<=1e-6)&&-fv>best,best=-fv;end; end
        catch; end
    end
    gam(k)=best;
end
figure('Position',[100 100 460 360]); b=bar(gam,0.5,'FaceColor',[0.2 0.45 0.75]);
hold on; yline(1,'--r','LineWidth',1.5);
text(2.2,1.03,'spec met (\gamma\geq1)','Color','r');
set(gca,'XTickLabel',mats(:,1)); ylabel('\gamma  (spec satisfaction ratio)');
title('Material comparison (d_c\geq0.5mm)'); ylim([0 1.2]); grid on;
for k=1:3, text(k,gam(k)+0.03,sprintf('%.3f',gam(k)),'HorizontalAlignment','center'); end
exportgraphics(gcf,fullfile(FIGDIR,'fig1_material_gamma.png'),'Resolution',200);

%% ===== Fig 2. d_c 민감도 (Al, surrogate) =====
dc=[0.50 0.40 0.35 0.30 0.26 0.22]; gdc=[0.753 0.874 0.955 1.053 1.151 1.273];
figure('Position',[100 100 480 360]);
plot(dc,gdc,'-o','LineWidth',1.8,'MarkerFaceColor','b'); hold on;
yline(1,'--r','LineWidth',1.5);
xline(0.32,':k','LineWidth',1.2);
text(0.33,0.6,'threshold \approx0.32mm','Rotation',90);
fill([0.2 0.32 0.32 0.2],[0 0 1.4 1.4],[1 0.9 0.9],'EdgeColor','none','FaceAlpha',0.4);
plot(dc,gdc,'-o','LineWidth',1.8,'MarkerFaceColor','b');
xlabel('coil wire diameter d_c [mm]'); ylabel('\gamma');
title('Spec satisfaction vs wire diameter (Al6061)'); grid on; xlim([0.2 0.5]); ylim([0.6 1.35]);
text(0.43,0.78,'d_c\geq0.5: not met','Color',[0.5 0 0]);
exportgraphics(gcf,fullfile(FIGDIR,'fig2_dc_sensitivity.png'),'Resolution',200);

%% ===== Fig 3. B_g surrogate vs FEMM (parity) =====
T=readtable(fullfile('..','femm','bg_grid.csv'));
Bfem=T.Bg_mT; Bsur=zeros(height(T),1);
for i=1:height(T)
    rm=T.r_m(i);hy=T.h_y1(i);tg=T.t_g(i);
    Bsur(i)=Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg];
end
R2=1-sum((Bsur-Bfem).^2)/sum((Bfem-mean(Bfem)).^2);
figure('Position',[100 100 400 380]);
plot(Bfem,Bsur,'o','MarkerFaceColor',[0.2 0.6 0.3],'MarkerEdgeColor','k'); hold on;
lim=[min(Bfem)-20 max(Bfem)+20]; plot(lim,lim,'--k','LineWidth',1.2);
xlabel('FEMM B_g [mT]'); ylabel('Surrogate B_g [mT]');
title(sprintf('B_g surrogate fit (R^2=%.3f)',R2)); axis equal; xlim(lim); ylim(lim); grid on;
exportgraphics(gcf,fullfile(FIGDIR,'fig3_bg_parity.png'),'Resolution',200);

%% ===== Fig 4. B_g vs t_g (설계 레버, opt 형상) =====
tgs=linspace(1,3,50); rm=6.23; hy=5.1;
Bsweep=arrayfun(@(tg)Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg],tgs);
tg_f=[1.5 2.0 2.5 3.0]; Bg_f=[382.2 338.8 306.2 276.9];   % FEMM 실측
figure('Position',[100 100 480 360]);
plot(tgs,Bsweep,'-','LineWidth',1.8,'Color',[0.2 0.45 0.75]); hold on;
plot(tg_f,Bg_f,'rs','MarkerFaceColor','r','MarkerSize',8);
xlabel('air-gap t_g [mm]'); ylabel('B_g [mT]'); grid on;
legend('surrogate','FEMM','Location','northeast');
title('Smaller gap \rightarrow higher B_g (design lever)');
exportgraphics(gcf,fullfile(FIGDIR,'fig4_bg_tg_lever.png'),'Resolution',200);

%% ===== Fig 5. stroke-f1 trade-off (binding constraint) =====
Kf=3.392; meff=0.0147;                 % 권장설계 근방(고정 근사)
kx=linspace(2000,30000,200);
stroke=Kf*2./kx*1e3;                   % mm @2A
f1=1/(2*pi)*sqrt(kx/meff);
figure('Position',[100 100 480 360]);
plot(stroke,f1,'-','LineWidth',1.8,'Color',[0.2 0.45 0.75]); hold on;
xline(1,'--r'); yline(100,'--r');
fill([1 3 3 1],[100 100 300 300],[0.85 1 0.85],'EdgeColor','none','FaceAlpha',0.4);
text(1.5,250,'spec region','Color',[0 0.5 0]);
plot(1.05,105,'kp','MarkerFaceColor','y','MarkerSize',14);
text(1.08,118,'d_c=0.3 design');
xlabel('max stroke @2A [mm]'); ylabel('1st resonance f_1 [Hz]');
title('Stroke vs resonance trade-off'); grid on; xlim([0 3]); ylim([40 300]);
exportgraphics(gcf,fullfile(FIGDIR,'fig5_tradeoff.png'),'Resolution',200);

%% ===== Fig 6. 최종설계 스펙 충족 (d_c=0.3) =====
figure('Position',[100 100 520 360]);
norm_req=[1 1 1];
X=[1.053/1, 105.3/100, 1.90/2];  % achieved/limit (stroke,f1: >=1 good; current: <=1 good)
bar([norm_req; X]','grouped');
set(gca,'XTickLabel',{'stroke/1mm','f_1/100Hz','I/2A'});
legend('requirement','achieved','Location','northwest');
ylabel('normalized'); title('Final design (d_c=0.3) vs spec'); grid on; ylim([0 1.3]);
yline(1,'--k');
exportgraphics(gcf,fullfile(FIGDIR,'fig6_final_spec.png'),'Resolution',200);

fprintf('생성 완료: fig1~fig6 (%s)\n', FIGDIR);
fprintf('material gamma: Al=%.3f, steel=%.3f, BeCu=%.3f\n', gam);

%% ===== 로컬 함수 (surrogate 기반 perf/nlcon) =====
function r=perf(xv,mat,Cs)
    t=xv(1)*1e-3;L=xv(2)*1e-3;b=xv(3)*1e-3;r_m=xv(4)*1e-3;h_m=xv(5)*1e-3;h_y1=xv(6)*1e-3;t_g=xv(7)*1e-3;d_c=xv(8)*1e-3;
    c_c=0.3e-3;t_b=0.3e-3;r_p=1.15;n_sym=2;I_max=2; a_s=15e-3;d_b=40e-3;h_i=5e-3; t_y3=1e-3;
    w_c=t_g-2*c_c;h_c=h_y1;
    rm=r_m*1e3;hy=h_y1*1e3;tg=t_g*1e3;
    Bg=(Cs*[1;rm;hy;tg;rm^2;hy^2;tg^2;rm*hy;rm*tg;hy*tg])/1000;
    A_g=2*pi*(r_m+t_g/2)*h_y1; By1=Bg*A_g/(pi*r_m^2);
    dcs=r_p*d_c; n_half=(h_c/dcs)*((2/sqrt(3))*(max(w_c,0)/dcs)); L_m=2*pi*(r_m+c_c+t_b+w_c/2);
    n_eff=min(n_sym*(h_y1+0.08*h_m)/h_c*n_half, n_sym*n_half); Kf=n_eff*Bg*L_m;
    m_coil=8960*(n_sym*n_half)*(pi/4)*d_c^2*L_m;
    Ia=b*t^3/12; k_x=24*mat.E*Ia/L^3;
    m_eff=mat.rho*a_s^2*b+m_coil+0.25*2*mat.rho*d_b*h_i*b+(1/3)*8*mat.rho*L*t*b;
    f1=(1/(2*pi))*sqrt(k_x/m_eff); x_max=Kf*I_max/k_x; sigma=3*mat.E*t*(x_max/2)/L^2;
    r.f1=f1;r.x_max=x_max;r.sigma=sigma;r.By1=By1;r.Bg=Bg;r.w_c=w_c;r.dcs=dcs;
    r.VCMd=2*(r_m+t_g+t_y3); r.footY=a_s+2*(2*L+h_i);
end
function [c,ceq]=nlc(xv,mat,Cs)
    r=perf(xv,mat,Cs); g=xv(9);
    c=[g*1e-3-r.x_max; g*100-r.f1; r.sigma-0.3*mat.Sy; r.By1-1.8; r.VCMd-20e-3; r.footY-200e-3; r.dcs-r.w_c; r.Bg-0.6];
    ceq=[];
end
