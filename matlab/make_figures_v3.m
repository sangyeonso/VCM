%% make_figures_v3.m - 전류 상향 최종안 figure (신규 current 스윕 + 최종 KPI 갱신)
clear; clc; close all;
fd=fullfile('..','figure'); blue=[0.15 0.39 0.92]; red=[0.80 0.20 0.20]; grn=[0.20 0.55 0.30]; gry=[0.55 0.55 0.6];

%% fig_v2_current : 전류 I -> γ (좌) + J (우)
I=[2.0 2.3 2.5 2.7 3.0]; G=[0.926 0.970 0.997 1.023 1.060]; J=[10.2 11.7 12.7 13.8 15.3];
f=figure('Position',[100 100 560 380],'Color','w');
yyaxis left; plot(I,G,'-o','Color',blue,'LineWidth',2,'MarkerFaceColor',blue); ylabel('\gamma (충족률)'); ylim([0.9 1.1]);
yline(1.0,'--','Color',gry,'LineWidth',1.2);
hold on; plot(2.7,1.023,'p','MarkerSize',16,'MarkerFaceColor',grn,'MarkerEdgeColor','k');
yyaxis right; plot(I,J,'-s','Color',red,'LineWidth',1.6,'MarkerFaceColor',red); ylabel('전류밀도 J [A/mm^2]'); ylim([9 16]);
yline(10,':','Color',red,'LineWidth',1.2);
set(gca,'FontSize',11); xlabel('구동 전류 I [A]  (원 스펙 2A)');
title({'전류 상향 → 단일 VCM으로 γ\geq1','I\approx2.7A에서 전 KPI 충족 (J=13.8, 총 P=1.1W)'},'FontSize',10);
text(2.72,1.023,'  채택 (2.7A, γ=1.02)','Color',grn,'FontSize',9,'FontWeight','bold');
exportgraphics(f,fullfile(fd,'fig_v2_current.png'),'Resolution',200);

%% fig_v2_final : 최종 KPI (2.7A 단일VCM) — 둘 다 충족
f=figure('Position',[100 100 560 380],'Color','w');
ach=[1.023 102.3]; tgt=[1.0 100]; rel=ach./tgt;
hb=bar([rel;[1 1]]','grouped'); hb(1).FaceColor=grn; hb(2).FaceColor=gry;
set(gca,'XTickLabel',{'행정 (1.02/1.0 mm)','공진 (102/100 Hz)'},'FontSize',10);
ylabel('목표 대비 비율'); ylim([0 1.2]); yline(1.0,'--k','목표','LineWidth',1.2);
legend({'달성 (\gamma=1.02)','목표'},'Location','south','FontSize',9);
title({'최종 설계: 단일 VCM, I=2.7A, DCP+Al7075 (d_c=0.5)','\gamma=1.02 — 전 성능 KPI 충족 (전류 2A\rightarrow2.7A 완화)'},'FontSize',9.5);
exportgraphics(f,fullfile(fd,'fig_v2_final.png'),'Resolution',200);

fprintf('생성: fig_v2_current.png, fig_v2_final.png (갱신)\n');
