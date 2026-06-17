"""
vcm_femm_halbach.py - 고자속 토폴로지 B_g 배율(beta) FEMM 실측
  baseline: 축방향 자석 + 철 요크 (현재 모델)
  radial  : 방사형 자석 링 (Halbach식 직접 radial gap drive)
  같은 Ø<=20 envelope에서 B_g 비교 -> beta = B_radial / B_baseline
"""
import os, math, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
import femm
TMP=r"C:\femm_vcm\hb.fem"; os.makedirs(r"C:\femm_vcm",exist_ok=True)
MU0=4*math.pi*1e-7; BR=1.43; MUY=1400

def setup(name):
    femm.newdocument(0); femm.mi_probdef(0,"millimeters","axi",1e-8,0,30)
    femm.mi_getmaterial("Air"); Hc=BR/(MU0*1.05)
    femm.mi_addmaterial("Nd",1.05,1.05,Hc,0,0,0,0,0,0,0,0)
    femm.mi_addmaterial("Fe",MUY,MUY,0,0,0,0,0,0,0,0,0)
def rect(r1,z1,r2,z2):
    femm.mi_addnode(r1,z1);femm.mi_addnode(r2,z1);femm.mi_addnode(r2,z2);femm.mi_addnode(r1,z2)
    femm.mi_addsegment(r1,z1,r2,z1);femm.mi_addsegment(r2,z1,r2,z2)
    femm.mi_addsegment(r2,z2,r1,z2);femm.mi_addsegment(r1,z2,r1,z1)
def blk(r,z,mat,md,g):
    femm.mi_addblocklabel(r,z);femm.mi_selectlabel(r,z)
    femm.mi_setblockprop(mat,1,0,"<None>",md,g,0);femm.mi_clearselected()
def gapB(rmid,z1,z2):
    femm.mi_zoomnatural();femm.mi_saveas(TMP);femm.mi_analyze(1);femm.mi_loadsolution()
    N=21;s=0.0
    for i in range(N):
        zz=z1+(z2-z1)*(i+0.5)/N; s+=abs(femm.mo_getpointvalues(rmid,zz)[1])
    return s/N

def baseline(r_m=6.0,h_m=8.0,h_y1=4.0,t_g=1.8,t_y3=0.8,h_y2=4.0):
    setup("base"); z0,z1=0.0,h_y2; z2=z1+h_m; z3=z2+h_y1
    r_si=r_m+t_g; r_so=r_m+t_g+t_y3
    rect(0,z0,r_so,z1);rect(0,z1,r_m,z2);rect(0,z2,r_m,z3);rect(r_si,z1,r_so,z3)
    Ro=4*r_so; femm.mi_makeABC(7,Ro,0,(z0+z3)/2,0)
    blk(r_m*0.5,(z0+z1)/2,"Fe",0,1);blk(r_m*0.5,(z1+z2)/2,"Nd",90,2)
    blk(r_m*0.5,(z2+z3)/2,"Fe",0,1);blk((r_si+r_so)/2,(z1+z3)/2,"Fe",0,1)
    blk(Ro*0.8,(z0+z3)/2,"Air",0,0)
    return gapB(r_m+t_g/2, z2, z3)

def radial(r_core=3.0,r_m=6.0,t_g=1.8,t_y=1.2,h=8.0,tp=3.0):
    # 방사형 자석: core(Fe) | magnet(radial +r) | gap | outer pole(Fe), 상하 plate(Fe)로 폐로
    setup("rad"); z1=tp; z2=tp+h; H=z2+tp
    r_si=r_m+t_g; r_so=r_si+t_y
    # 상/하 plate (core~outer pole 연결)
    rect(0,0,r_so,z1); rect(0,z2,r_so,H)
    rect(0,z1,r_core,z2)          # core
    rect(r_core,z1,r_m,z2)        # radial magnet
    rect(r_si,z1,r_so,z2)         # outer pole
    Ro=4*r_so; femm.mi_makeABC(7,Ro,0,H/2,0)
    blk(r_so*0.5,z1*0.5,"Fe",0,1); blk(r_so*0.5,(z2+H)/2,"Fe",0,1)   # plates
    blk(r_core*0.5,(z1+z2)/2,"Fe",0,1)                                # core
    blk((r_core+r_m)/2,(z1+z2)/2,"Nd",0,2)                            # radial magnet (+r)
    blk((r_si+r_so)/2,(z1+z2)/2,"Fe",0,1)                             # outer pole
    blk((r_m+r_si)/2,(z1+z2)/2,"Air",0,0)                             # gap (코일 영역)
    blk(Ro*0.8,H/2,"Air",0,0)
    return gapB(r_m+t_g/2, z1, z2)

def main():
    femm.openfemm(1)
    Bb=baseline()
    print(f"baseline(축+철)  B_g = {Bb*1000:.0f} mT  (필요 beta=1.6)")
    print(f"{'r_core':>7}{'h_mag':>7}{'Bg[mT]':>9}{'beta':>7}")
    best=0
    for rc in [1.0,1.5,2.0,3.0]:
        for h in [8,12,16]:
            try:
                Br=radial(r_core=rc,h=h); be=Br/Bb
                if be>best: best=be
                print(f"{rc:>7.1f}{h:>7.0f}{Br*1000:>9.0f}{be:>7.2f}")
            except Exception as e:
                print(f"{rc:>7.1f}{h:>7.0f}   실패 {str(e)[:20]}")
    femm.closefemm()
    print(f"\n최대 beta(방사자석) = {best:.2f}  vs  필요 1.6")

if __name__=="__main__": main()
