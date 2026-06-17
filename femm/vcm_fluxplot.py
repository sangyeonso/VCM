"""
vcm_fluxplot.py  -  권장 최종설계 VCM의 FEMM 자속밀도 분포 그림 저장 (보고서용)
실행: (GPSM)  python vcm_fluxplot.py   -> figure/fig8_femm_flux.png
주의: FEMM 창이 보이게 열려야 savebitmap이 화면을 캡처함(openfemm(0)).
"""
import os, math, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
import femm
from PIL import Image

TMP=r"C:\femm_vcm\flux_tmp.fem"; os.makedirs(r"C:\femm_vcm",exist_ok=True)
HERE=os.path.dirname(__file__); FIG=os.path.join(HERE,"..","figure")
BMP=r"C:\femm_vcm\flux.bmp"; PNG=os.path.join(FIG,"fig8_femm_flux.png")
MU0=4*math.pi*1e-7; BR=1.17; MU_YOKE=1400

# 권장 최종설계(d_c=0.3) VCM 형상 [mm]
r_m,h_m,h_y1,t_g,t_y3,h_y2 = 6.5,25.0,5.8,2.5,1.0,3.0
z0,z1=0.0,h_y2; z2=z1+h_m; z3=z2+h_y1; r_si=r_m+t_g; r_so=r_m+t_g+t_y3

femm.openfemm(0); femm.main_resize(480,1000)   # 세로형 창(VCM이 길쭉)
femm.newdocument(0)
femm.mi_probdef(0,"millimeters","axi",1e-8,0,30)
femm.mi_getmaterial("Air")
femm.mi_addmaterial("NdFeB117",1.05,1.05,BR/(MU0*1.05),0,0,0,0,0,0,0,0)
femm.mi_addmaterial("Yoke",MU_YOKE,MU_YOKE,0,0,0,0,0,0,0,0,0)
def rect(r1,a,r2,b):
    femm.mi_addnode(r1,a);femm.mi_addnode(r2,a);femm.mi_addnode(r2,b);femm.mi_addnode(r1,b)
    femm.mi_addsegment(r1,a,r2,a);femm.mi_addsegment(r2,a,r2,b)
    femm.mi_addsegment(r2,b,r1,b);femm.mi_addsegment(r1,b,r1,a)
rect(0,z0,r_so,z1); rect(0,z1,r_m,z2); rect(0,z2,r_m,z3); rect(r_si,z1,r_so,z3)
for rr,zz in [(r_m,z2),(r_si,z2),(r_m,z3),(r_si,z3)]: femm.mi_addnode(rr,zz)
femm.mi_addsegment(r_m,z2,r_si,z2); femm.mi_addsegment(r_m,z3,r_si,z3)
Rout=2.5*max(r_so,z3); femm.mi_makeABC(7,Rout,0,(z0+z3)/2,0)
def blk(r,z,mat,md,g):
    femm.mi_addblocklabel(r,z);femm.mi_selectlabel(r,z)
    femm.mi_setblockprop(mat,1,0,"<None>",md,g,0);femm.mi_clearselected()
blk(r_m*0.5,(z0+z1)/2,"Yoke",0,1); blk(r_m*0.5,(z1+z2)/2,"NdFeB117",90,2)
blk(r_m*0.5,(z2+z3)/2,"Yoke",0,1); blk((r_si+r_so)/2,(z1+z3)/2,"Yoke",0,1)
blk((r_m+r_si)/2,(z2+z3)/2,"Air",0,0); blk((r_m+r_si)/2,(z1+z2)/2,"Air",0,0)
blk(Rout*0.8,(z0+z3)/2,"Air",0,0)

femm.mi_zoomnatural(); femm.mi_saveas(TMP); femm.mi_analyze(1); femm.mi_loadsolution()
# 자속밀도 컬러맵 + 자속선
femm.mo_hidegrid()
femm.mo_showdensityplot(1,0,1.2,0,"bmag")   # legend, color, upper=1.2T, lower=0, |B|
femm.mo_zoom(-1,-2,13,z3+3)                 # VCM 영역 확대
femm.mo_savebitmap(BMP)
femm.closefemm()

# BMP -> PNG 변환(+여백 약간 트림)
img=Image.open(BMP).convert("RGB")
img.save(PNG)
print("저장:",PNG,"  크기:",img.size)
