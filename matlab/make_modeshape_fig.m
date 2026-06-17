%% make_modeshape_fig.m - 1차 모드형상 (올바른 DCP 구조, 정정본)
%  구조: 중앙 플랫폼 + 전폭 상하 바 + 양옆 4 지점핀 + 16 힌지 (fig0와 동일)
%  1차 모드 = X 병진: 플랫폼 +D, 중간바 +D/2, 지점핀 고정, blade는 guided S자 변형
clear; clc; close all;
figure('Position',[100 100 560 560]); hold on; axis equal off;
gry=[0.62 0.62 0.66]; blue=[0.10 0.30 0.90];
D=9;  S=@(s) 3*s.^2-2*s.^3;            % guided-beam S자 형상, D=모드 스케일(과장)

% blade 정의: {x0, y_bot, y_top, shift_bot, shift_top}  (shift: 핀0 / 바 D/2 / 플랫폼 D)
B={ {46,60,86,D,D/2},{64,60,86,D,D/2}, {46,14,40,D/2,D},{64,14,40,D/2,D}, ...   % 내부
    {28,60,86,0,D/2},{82,60,86,0,D/2}, {28,14,40,D/2,0},{82,14,40,D/2,0} };     % 외부

% ===== 변형 전(undeformed, 회색) =====
rectangle('Position',[24 86 62 5],'EdgeColor',gry,'FaceColor','none');
rectangle('Position',[24  9 62 5],'EdgeColor',gry,'FaceColor','none');
rectangle('Position',[40 40 30 20],'Curvature',0.2,'EdgeColor',gry,'FaceColor','none');
for i=1:8, b=B{i}; plot([b{1} b{1}],[b{2} b{3}],'-','Color',gry,'LineWidth',1); end

% ===== 1차 모드(deformed, 파랑) =====
rectangle('Position',[24+D/2 86 62 5],'EdgeColor',blue,'FaceColor','none','LineWidth',1.5);
rectangle('Position',[24+D/2  9 62 5],'EdgeColor',blue,'FaceColor','none','LineWidth',1.5);
rectangle('Position',[40+D 40 30 20],'Curvature',0.2,'EdgeColor',blue,'FaceColor','none','LineWidth',1.8);
for i=1:8
  b=B{i}; x0=b{1}; yb=b{2}; yt=b{3}; sb=b{4}; st=b{5};
  s=linspace(0,1,30); y=yb+(yt-yb)*s; x=x0+sb+(st-sb)*S(s);
  plot(x,y,'-','Color',blue,'LineWidth',2.2);
end
% 지점핀(고정) 표시
for p=[28 60;82 60;28 40;82 40]', plot(p(1),p(2),'^','Color',gry,'MarkerFaceColor',gry,'MarkerSize',6); end

text(55,98,'1st mode = X-translation (DCP)','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
text(55,94,'f_1 \approx 93 Hz (최종 설계, 컴플라이언스 행렬법)','HorizontalAlignment','center','FontSize',9,'Color',[0.3 0.35 0.4]);
% 범례
plot([20 30],[3 3],'-','Color',gry,'LineWidth',1); text(32,3,'undeformed','FontSize',9,'Color',gry,'VerticalAlignment','middle');
plot([62 72],[3 3],'-','Color',blue,'LineWidth',2.2); text(74,3,'1st mode (scaled)','FontSize',9,'Color',blue,'VerticalAlignment','middle');

xlim([4 114]); ylim([0 100]);
exportgraphics(gcf,fullfile('..','figure','fig7_fea_modeshape.png'),'Resolution',200);
fprintf('저장: figure/fig7_fea_modeshape.png (올바른 DCP 구조, f1=93Hz)\n');
