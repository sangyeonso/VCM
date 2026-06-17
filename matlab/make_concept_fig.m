%% make_concept_fig.m - 개념도(DCP 유연기구 + VCM) PNG 생성 (LaTeX 삽입용)
clear; clc; close all;
figure('Position',[100 100 640 540]); hold on; axis equal off;
blue=[0.15 0.39 0.92]; fixc=[0.90 0.92 0.95]; soft=[0.86 0.91 0.99]; edg=[0.5 0.55 0.6];

% 고정 베이스(ground) 상·하
rectangle('Position',[8 6 78 6],'FaceColor',fixc,'EdgeColor',edg);
rectangle('Position',[8 88 78 6],'FaceColor',fixc,'EdgeColor',edg);
text(47,99,'Fixed base (Ground)','HorizontalAlignment','center','FontSize',9,'Color',[0.4 0.45 0.5]);

% 중간 스테이지
rectangle('Position',[16 28 58 5],'FaceColor',[0.93 0.95 0.98],'EdgeColor',[0.58 0.64 0.72]);
rectangle('Position',[16 67 58 5],'FaceColor',[0.93 0.95 0.98],'EdgeColor',[0.58 0.64 0.72]);
text(45,30.5,'Intermediate','HorizontalAlignment','center','FontSize',8);
text(45,69.5,'Intermediate','HorizontalAlignment','center','FontSize',8);

% 무빙 스테이지
rectangle('Position',[12 45 64 10],'Curvature',0.12,'FaceColor',soft,'EdgeColor',blue,'LineWidth',1.5);
text(36,50,'Moving Stage (coil)','HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.1 0.2 0.45]);

% blade 8개 (x=24, 62 ; 4단)
for x=[24 62]
  plot([x x],[12 28],'-','Color',blue,'LineWidth',2.5);
  plot([x x],[33 45],'-','Color',blue,'LineWidth',2.5);
  plot([x x],[55 67],'-','Color',blue,'LineWidth',2.5);
  plot([x x],[72 88],'-','Color',blue,'LineWidth',2.5);
end
text(15,40,'leaf blade \times8','FontSize',8,'Color',[0.4 0.45 0.5]);

% 치수 L, t (좌상단 blade)
plot([20 20],[72 88],'-','Color',[0.5 0.5 0.5]);
plot([19 21],[72 72],'-','Color',[0.5 0.5 0.5]); plot([19 21],[88 88],'-','Color',[0.5 0.5 0.5]);
text(17.5,80,'L','HorizontalAlignment','right','FontSize',10);
text(24,90.5,'t','HorizontalAlignment','center','FontSize',10);

% VCM: 코일(가동) + 자석/요크(고정)
rectangle('Position',[78 46 7 8],'FaceColor',blue,'EdgeColor',blue);
text(81.5,42.5,'coil','HorizontalAlignment','center','FontSize',8);
rectangle('Position',[88 43 15 14],'FaceColor',fixc,'EdgeColor',edg);
rectangle('Position',[90 50 10 4.5],'FaceColor',[0.93 0.55 0.55],'EdgeColor','none');
rectangle('Position',[90 45.5 10 4.5],'FaceColor',[0.58 0.68 0.95],'EdgeColor','none');
text(95,52.2,'N','HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
text(95,47.7,'S','HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
text(95.5,39.5,'magnet+yoke (fixed)','HorizontalAlignment','center','FontSize',7.5,'Color',[0.4 0.45 0.5]);
text(95.5,60,'VCM','HorizontalAlignment','center','FontSize',8.5,'Color',[0.3 0.35 0.4],'FontWeight','bold');

% X 운동 화살표
plot([30 46],[38 38],'-','Color',blue,'LineWidth',2);
plot(30,38,'<','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',8);
plot(46,38,'>','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',8);
text(38,34.5,'X  (1-DOF, >1mm)','HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.1 0.2 0.45]);

% 좌표축
plot([100 108],[12 12],'-k','LineWidth',1); plot(108,12,'>k','MarkerFaceColor','k','MarkerSize',5);
plot([100 100],[12 20],'-k','LineWidth',1); plot(100,20,'^k','MarkerFaceColor','k','MarkerSize',5);
text(109,12,'X','FontSize',8); text(100,22,'Y','HorizontalAlignment','center','FontSize',8);

xlim([4 114]); ylim([2 104]);
exportgraphics(gcf,fullfile('..','figure','fig0_concept.png'),'Resolution',200);
fprintf('저장: figure/fig0_concept.png\n');
