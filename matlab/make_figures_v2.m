%% make_figures_v2.m - 정정된 분석 기반 figure 재생성 (전수감사 결과 반영)
%  계산 완료된 값들을 플롯. 저장: ../figure/fig_v2_*.png
clear; clc; close all;
fd = fullfile('..','figure');
blue=[0.15 0.39 0.92]; red=[0.80 0.20 0.20]; grn=[0.20 0.55 0.30]; gry=[0.55 0.55 0.6];

%% (1) 유연기구 수식 정정: 단순식 vs 검증된 행렬법  [fig_v2_method]
f=figure('Position',[100 100 560 380],'Color','w');
G=[0.873 0.947 0.905;   % Al7075: 단순24EI/L3, 행렬-단일, 행렬-DCP
   0.923 1.020 0.985];  % Mg
hb=bar(G,'grouped'); hb(1).FaceColor=gry; hb(2).FaceColor=blue; hb(3).FaceColor=grn;
set(gca,'XTickLabel',{'Al7075','Mg'},'FontSize',11); ylabel('\gamma (스펙 충족률)'); ylim([0.8 1.1]);
yline(1.0,'--k','\gamma=1 (충족)','LineWidth',1.2,'FontSize',10);
legend({'단순식 24EI/L^3 (틀림)','행렬법 단일','행렬법 DCP'},'Location','southoutside','Orientation','horizontal','FontSize',9);
title('유연기구 수식 정정 효과 (d_c=0.5, 발열·피로 포함)','FontSize',11);
exportgraphics(f,fullfile(fd,'fig_v2_method.png'),'Resolution',200);

%% (2) 발열 한계: J 제약이 d_c를 스펙(0.5)으로 밀고 γ를 가둠  [fig_v2_thermal]
f=figure('Position',[100 100 560 380],'Color','w');
Jm=[64 10.2 7]; dcopt=[0.20 0.505 0.603]; gam=[1.30 0.868 0.774];
yyaxis left;  plot(Jm,gam,'-o','Color',blue,'LineWidth',2,'MarkerFaceColor',blue); ylabel('\gamma'); ylim([0.6 1.4]);
yline(1.0,'--','Color',gry);
yyaxis right; plot(Jm,dcopt,'-s','Color',red,'LineWidth',2,'MarkerFaceColor',red); ylabel('최적 d_c [mm]'); ylim([0 0.7]);
yline(0.5,':','Color',red,'LineWidth',1.2);
set(gca,'XDir','reverse','FontSize',11); xlabel('허용 전류밀도 J_{max} [A/mm^2]  (왼쪽=발열무시, 오른쪽=현실)');
title({'발열(J) 한계가 d_c>0.5 스펙을 강제','d_c 무시 시 0.2mm로 몰리나(비현실), J\leq10이면 d_c\approx0.5'},'FontSize',10);
text(10.2,0.55,'스펙 d_c=0.5','Color',red,'FontSize',9);
exportgraphics(f,fullfile(fd,'fig_v2_thermal.png'),'Resolution',200);

%% (3) 재료 공동최적 랭킹 (DCP, d_c=0.5)  [fig_v2_material]
f=figure('Position',[100 100 560 380],'Color','w');
mats={'Mg','Al6061','Al7075','Ti-6Al-4V','BeCu','17-4PH','Spring','Maraging'};
gM=[0.985 0.912 0.905 0.827 0.744 0.752 0.751 0.746];
[gs,ix]=sort(gM);
hb=barh(gs); hb.FaceColor='flat';
for i=1:8, if gs(i)>=1, hb.CData(i,:)=grn; elseif gs(i)>=0.9, hb.CData(i,:)=blue; else, hb.CData(i,:)=gry; end; end
set(gca,'YTickLabel',mats(ix),'FontSize',10); xlabel('\gamma'); xlim([0.7 1.05]);
xline(1.0,'--k','\gamma=1','LineWidth',1.2);
title('재료 공동최적 랭킹 (DCP, d_c=0.5, 발열·피로 포함)','FontSize',11);
exportgraphics(f,fullfile(fd,'fig_v2_material.png'),'Resolution',200);

%% (4) 최종 설계 KPI (DCP + Al7075, d_c=0.5)  [fig_v2_final]
f=figure('Position',[100 100 560 380],'Color','w');
ach=[0.926 92.6]; tgt=[1.0 100];
rel=ach./tgt;   % 목표 대비
hb=bar([rel; [1 1]]','grouped'); hb(1).FaceColor=blue; hb(2).FaceColor=gry;
set(gca,'XTickLabel',{'행정 (0.93/1.0 mm)','공진 (92.6/100 Hz)'},'FontSize',10);
ylabel('목표 대비 비율'); ylim([0 1.2]); yline(1.0,'--k','목표','LineWidth',1.2);
legend({'달성 (\gamma=0.926)','목표'},'Location','south','FontSize',9);
title({'최종 설계 DCP+Al7075 (d_c=0.5 스펙준수)','\gamma=0.926 — 응력한계로 ~7% 미달 (경계 위)'},'FontSize',10);
exportgraphics(f,fullfile(fd,'fig_v2_final.png'),'Resolution',200);

fprintf('생성: fig_v2_method / fig_v2_thermal / fig_v2_material / fig_v2_final .png\n');
