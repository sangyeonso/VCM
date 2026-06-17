%% animate_flexure.m - 올바른 DCP 구조 X-스캐닝 애니메이션 (GIF)
%  구조: 중앙 플랫폼 M + 전폭 상하 바 + 양옆 4 지점핀 + 16 힌지 (검증 구조)
%  모션: 플랫폼 +d, 중간바 +d/2(compound), 링크 틸트, 지점핀 고정. (변위 과장)
clear; clc; close all;
gif=fullfile('..','figure','fig13_flexure_anim.gif');
blue=[0.15 0.39 0.92]; fixc=[0.90 0.92 0.95]; soft=[0.86 0.91 0.99]; edg=[0.5 0.55 0.6]; gry=[0.45 0.45 0.5];
A=7;                         % 화면상 변위 진폭(과장), 실제 ±1mm 대응
N=36; fig=figure('Position',[100 100 600 520],'Color','w');
% 고정 지점핀 좌표
pins=[28 60 -1; 82 60 1; 28 40 -1; 82 40 1];
for f=1:N
  d=A*sin(2*pi*(f-1)/N); di=d/2;          % 플랫폼 d, 중간바 di
  clf; hold on; axis equal off; xlim([4 114]); ylim([2 100]);
  % 전폭 상/하 중간바 (+di)
  rectangle('Position',[24+di 86 62 5],'FaceColor',[0.93 0.95 0.98],'EdgeColor',[0.58 0.64 0.72]);
  rectangle('Position',[24+di  9 62 5],'FaceColor',[0.93 0.95 0.98],'EdgeColor',[0.58 0.64 0.72]);
  % 중앙 플랫폼 M (+d)
  rectangle('Position',[40+d 40 30 20],'Curvature',0.25,'FaceColor',soft,'EdgeColor',blue,'LineWidth',1.6);
  text(55+d,50,'M','HorizontalAlignment','center','FontSize',15,'FontWeight','bold','Color',[0.1 0.2 0.45],'FontAngle','italic');
  % 링크 8개 (지점핀 고정 / 바 +di / 플랫폼 +d)
  segs={ {46+di,86,46+d,60},{64+di,86,64+d,60},{46+d,40,46+di,14},{64+d,40,64+di,14}, ... % 내부
         {28+di,86,28,60},  {82+di,86,82,60},  {28,40,28+di,14},  {82,40,82+di,14} };     % 외부
  for q=1:8
    plot([segs{q}{1} segs{q}{3}],[segs{q}{2} segs{q}{4}],'-','Color',blue,'LineWidth',2.4);
  end
  % 16 힌지
  HX=[28+di 46+di 64+di 82+di  46+d 64+d 46+d 64+d  28+di 46+di 64+di 82+di  28 82 28 82];
  HY=[86 86 86 86  60 60 40 40  14 14 14 14  60 60 40 40];
  plot(HX,HY,'o','MarkerFaceColor','w','MarkerEdgeColor',blue,'MarkerSize',5,'LineWidth',1.2);
  % 지점핀 (해치 삼각, 고정)
  for i=1:4
    px=pins(i,1);py=pins(i,2);dd=pins(i,3);
    patch(px+[0 -7*dd -7*dd],py+[0 4 -4],fixc,'EdgeColor',edg);
  end
  % VCM 코일(+d, 가동) + 자석(고정)
  plot([70+d 82+d],[50 50],'-','Color',blue,'LineWidth',1.5);
  rectangle('Position',[82+d 45 6 10],'FaceColor',blue,'EdgeColor',blue);
  rectangle('Position',[93 42 14 16],'FaceColor',fixc,'EdgeColor',edg);
  rectangle('Position',[95 50 10 5],'FaceColor',[0.93 0.55 0.55],'EdgeColor','none');
  rectangle('Position',[95 45 10 5],'FaceColor',[0.58 0.68 0.95],'EdgeColor','none');
  text(100,52.5,'N','HorizontalAlignment','center','FontSize',7,'Color','w','FontWeight','bold');
  text(100,47.5,'S','HorizontalAlignment','center','FontSize',7,'Color','w','FontWeight','bold');
  % X 변위 화살표
  plot([55 55+d],[30 30],'-','Color',[0.8 0.2 0.2],'LineWidth',2);
  plot(55+d,30,'o','MarkerFaceColor',[0.8 0.2 0.2],'MarkerEdgeColor',[0.8 0.2 0.2],'MarkerSize',5);
  title(sprintf('VCM-driven DCP flexure: X-scanning  (x = %+.2f mm)', d/A*1.0),'FontSize',11);
  drawnow;
  fr=getframe(fig); [ind,cm]=rgb2ind(frame2im(fr),256);
  if f==1, imwrite(ind,cm,gif,'gif','LoopCount',inf,'DelayTime',0.06);
  else,    imwrite(ind,cm,gif,'gif','WriteMode','append','DelayTime',0.06); end
  if f==10, exportgraphics(fig,fullfile('..','figure','fig13_frame.png'),'Resolution',150); end
end
fprintf('애니메이션 저장: %s (%d 프레임, 올바른 DCP 구조)\n', gif, N);
