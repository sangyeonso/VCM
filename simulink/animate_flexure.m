%% animate_flexure.m - DCP 유연기구 X-스캐닝 모션 애니메이션 (GIF)
%  stage 변위 d(t)=A sin, 중간stage d/2, blade는 guided S자 변형. (시인성 위해 변위 과장)
clear; clc; close all;
gif=fullfile('..','figure','fig13_flexure_anim.gif');
blue=[0.15 0.39 0.92]; fixc=[0.88 0.90 0.93]; soft=[0.86 0.91 0.99]; edg=[0.5 0.55 0.6];
A=14;                       % 화면상 변위 진폭(과장), 실제 ±1mm 대응
S=@(s) 3*s.^2-2*s.^3;       % guided-beam S자 형상
blade=@(xb,yb,xt,yt) deal( xb+(xt-xb)*S(linspace(0,1,24)), yb+(yt-yb)*linspace(0,1,24) );

N=36; fig=figure('Position',[100 100 560 520],'Color','w');
for f=1:N
  d=A*sin(2*pi*(f-1)/N);  di=d/2;          % stage 변위 d, 중간 di
  clf; hold on; axis equal off; xlim([-10 124]); ylim([2 104]);
  % 고정 ground (상/하)
  rectangle('Position',[8 6 78 6],'FaceColor',fixc,'EdgeColor',edg);
  rectangle('Position',[8 88 78 6],'FaceColor',fixc,'EdgeColor',edg);
  % 중간 stage (di), 무빙 stage (d)
  rectangle('Position',[16+di 28 58 5],'FaceColor',[0.93 0.95 0.98],'EdgeColor',[0.58 0.64 0.72]);
  rectangle('Position',[16+di 67 58 5],'FaceColor',[0.93 0.95 0.98],'EdgeColor',[0.58 0.64 0.72]);
  rectangle('Position',[12+d 45 64 10],'Curvature',0.12,'FaceColor',soft,'EdgeColor',blue,'LineWidth',1.5);
  text(40+d,50,'Moving Stage','HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.1 0.2 0.45]);
  % blade 8개: 아래로부터  ground(0)->inter(di)->stage(d)->inter(di)->ground(0)
  for xc=[24 62]
    segs={ {xc,12, xc+di,28}, {xc+di,33, xc+d,45}, {xc+d,55, xc+di,67}, {xc+di,72, xc,88} };
    for q=1:4
      [bx,by]=blade(segs{q}{1},segs{q}{2},segs{q}{3},segs{q}{4});
      plot(bx,by,'-','Color',blue,'LineWidth',2.4);
    end
  end
  % VCM 코일(가동, stage와 함께) + 고정 자석
  rectangle('Position',[78+d 46 7 8],'FaceColor',blue,'EdgeColor',blue);
  rectangle('Position',[92 43 15 14],'FaceColor',fixc,'EdgeColor',edg);
  rectangle('Position',[94 50 11 4.5],'FaceColor',[0.93 0.55 0.55],'EdgeColor','none');
  rectangle('Position',[94 45.5 11 4.5],'FaceColor',[0.58 0.68 0.95],'EdgeColor','none');
  text(99,52.2,'N','HorizontalAlignment','center','FontSize',7,'Color','w','FontWeight','bold');
  text(99,47.7,'S','HorizontalAlignment','center','FontSize',7,'Color','w','FontWeight','bold');
  % X 화살표 + 변위 표시
  plot([54 54+d],[38 38],'-','Color',[0.8 0.2 0.2],'LineWidth',2);
  plot(54+d,38,'o','MarkerFaceColor',[0.8 0.2 0.2],'MarkerEdgeColor',[0.8 0.2 0.2],'MarkerSize',5);
  title(sprintf('VCM-driven DCP flexure: X-scanning  (x = %+.2f mm)', d/A*1.0),'FontSize',11);
  drawnow;
  % GIF 기록
  fr=getframe(fig); [ind,cm]=rgb2ind(frame2im(fr),256);
  if f==1, imwrite(ind,cm,gif,'gif','LoopCount',inf,'DelayTime',0.06);
  else,    imwrite(ind,cm,gif,'gif','WriteMode','append','DelayTime',0.06); end
end
fprintf('애니메이션 저장: %s (%d 프레임)\n', gif, N);
