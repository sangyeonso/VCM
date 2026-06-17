%% make_concept_fig.m - 개념도(올바른 DCP 구조 + VCM) PNG 생성
%  구조: 중앙 플랫폼 M + 전폭 상/하 중간바(m) + 양옆 4 지점핀 + 16 힌지 (강의/검증 구조)
clear; clc; close all;
figure('Position',[100 100 660 540]); hold on; axis equal off;
blue=[0.15 0.39 0.92]; fixc=[0.90 0.92 0.95]; soft=[0.86 0.91 0.99]; edg=[0.5 0.55 0.6]; gry=[0.45 0.45 0.5];

% 전폭 상/하 중간바 (m)
rectangle('Position',[24 86 62 5],'FaceColor',[0.93 0.95 0.98],'EdgeColor',[0.58 0.64 0.72]);
rectangle('Position',[24  9 62 5],'FaceColor',[0.93 0.95 0.98],'EdgeColor',[0.58 0.64 0.72]);
text(55,93.5,'중간 바 (intermediate, m)','HorizontalAlignment','center','FontSize',8,'Color',gry);
text(55,6.5 ,'중간 바 (intermediate, m)','HorizontalAlignment','center','FontSize',8,'Color',gry);

% 중앙 무빙 플랫폼 M
rectangle('Position',[40 40 30 20],'Curvature',0.25,'FaceColor',soft,'EdgeColor',blue,'LineWidth',1.6);
text(55,50,'M','HorizontalAlignment','center','FontSize',16,'FontWeight','bold','Color',[0.1 0.2 0.45],'FontAngle','italic');

% 링크 8개 (수직, blade)
L={[46 86 46 60],[64 86 64 60],[46 40 46 14],[64 40 64 14],...   % 내부(바<->플랫폼)
   [28 86 28 60],[82 86 82 60],[28 14 28 40],[82 14 82 40]};     % 외부(바<->지점핀)
for i=1:8, plot([L{i}(1) L{i}(3)],[L{i}(2) L{i}(4)],'-','Color',blue,'LineWidth',2.4); end

% 16 힌지
HX=[28 46 64 82  46 64 46 64  28 46 64 82  28 82 28 82];
HY=[86 86 86 86  60 60 40 40  14 14 14 14  60 60 40 40];
plot(HX,HY,'o','MarkerFaceColor','w','MarkerEdgeColor',blue,'MarkerSize',5,'LineWidth',1.2);

% 양옆 4 지점핀 (해치 삼각) : (28,60)(82,60)(28,40)(82,40)
pins=[28 60 -1; 82 60 1; 28 40 -1; 82 40 1];
for i=1:4
  px=pins(i,1);py=pins(i,2);d=pins(i,3);
  patch(px+[0 -7*d -7*d],py+[0 4 -4],fixc,'EdgeColor',edg);
  for k=-1:1, plot(px-7*d+[0 -2.5*d],py+4*k+[0 -2.5],'-','Color',gry,'LineWidth',0.7); end
end
text(13,50,'지점핀 \times4','HorizontalAlignment','center','FontSize',8,'Color',gry,'Rotation',90);

% L 치수 (좌상단 외부 링크)
plot([20 20],[60 86],'-','Color',gry); plot([19 21],[60 60],'-','Color',gry); plot([19 21],[86 86],'-','Color',gry);
text(17,73,'L','HorizontalAlignment','right','FontSize',11);

% VCM: 코일(가동, 플랫폼에 연결) + 자석/요크(고정)
plot([70 82],[50 50],'-','Color',blue,'LineWidth',1.5);            % 연결 로드
rectangle('Position',[82 45 6 10],'FaceColor',blue,'EdgeColor',blue);
text(85,42,'coil','HorizontalAlignment','center','FontSize',8,'Color',[0.1 0.2 0.45]);
rectangle('Position',[91 42 15 16],'FaceColor',fixc,'EdgeColor',edg);
rectangle('Position',[93 50 11 5],'FaceColor',[0.93 0.55 0.55],'EdgeColor','none');
rectangle('Position',[93 45 11 5],'FaceColor',[0.58 0.68 0.95],'EdgeColor','none');
text(98.5,52.5,'N','HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
text(98.5,47.5,'S','HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
text(98.5,60.5,'VCM','HorizontalAlignment','center','FontSize',9,'Color',[0.3 0.35 0.4],'FontWeight','bold');
text(98.5,39,'magnet+yoke (고정)','HorizontalAlignment','center','FontSize',7,'Color',gry);

% X 운동 화살표
plot([47 63],[30 30],'-','Color',blue,'LineWidth',2);
plot(47,30,'<','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',7);
plot(63,30,'>','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',7);
text(55,26,'x  (1-DOF, >1mm)','HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.1 0.2 0.45]);

text(55,70,'16 leaf hinges','HorizontalAlignment','center','FontSize',8,'Color',gry);

xlim([4 114]); ylim([2 100]);
exportgraphics(gcf,fullfile('..','figure','fig0_concept.png'),'Resolution',200);
fprintf('저장: figure/fig0_concept.png (올바른 DCP 구조)\n');
