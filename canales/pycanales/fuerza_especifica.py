import numpy as np
from math import log,atan,acos,pi,sin,cos 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon
################################################
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
################################
from matplotlib.patches import ConnectionPatch
from matplotlib.gridspec import GridSpec
#######################
def makemat(a,f,c):
    mat_resul=[]
    for i in range(f):
        mat_resul.append([])
        for j in range(c):
            mat_resul[i].append(a)
    return mat_resul



def E(J0,r,t):
    res=J0**4+((5*t+2)/(2))*J0**3+(((3*t+2)*(t+1))/(2))*J0**2+(t**2/2+(t-6*r)*(t+1))*J0-6*r*(t+1)**2
    return res

def dE(J0,r,t):
    res=4*J0**3+3*((5*t+2)/(2))*J0**2+2*(((3*t+2)*(t+1))/(2))*J0+(t**2/2+(t-6*r)*(t+1))
    return res
#####################

def tabulador_ej(a,b,n,r,t):
    inc=(b-a)/n
    x=[a+inc*i for i in range(n+1)]
    y=[E(a+inc*i,r,t) for i in range(n+1)]
    return [x,y]
def FET(Q,b,z1,z2,y):
    g=9.807
    A=0.5*(z1+z2)*y**2+b*y
    ym=(1/3+1/6*((b*y)/A))*y
    res=(Q**2)/(g*A)+ym*A
    return res
def EET(Q,b,z1,z2,alfa,y):
    g=9.807
    A=0.5*(z1+z2)*y**2+b*y
    V=Q/A
    res=y+alfa*(V**2)/(2*g)
    return res
################################################
class Resalto_r:
    def __init__(self,Q,y,b):
        self.Q=Q
        self.y=y
        self.b=b
    def results(self,decimal_number=6,tol=0.0001,show_results="TRUE"):
    ###########################
        A=self.b*self.y
        V=self.Q/A
        g=9.807
        T=self.b
        F=(V)/(g*(A/T))**0.5
        y2=(self.y/2)*((1+8*F**2)**0.5-1)
        if (show_results=="TRUE"):
            print("=================================================================")
            if (y2>self.y):
                print("Tirante supercrítico    (y1)[m] :",self.y)
                print("Tirante subcrítico      (y2)[m] :",y2)
            if (y2<self.y):
                print("Tirante supercrítico    (y1)[m] :",y2)
                print("Tirante subcrítico      (y2)[m] :",self.y)
        return y2        
    def show_FY(self,lis_py=[],lis_col=["blue","red","orange","gray"],posx=0.80,posy=0.5,tam=0.1,dimy=3*0.6):
        n=50
        fac=tam
        inc=(dimy)/n
        Y=[inc*(i+1) for i in range(n)]
        asintota=[]
        for i in range(n):
            A=self.b*Y[i]
            ym=0.5*Y[i]
            asintota.append(ym*A)    
        lis_F=[]
        ######################
        lis_Q=[self.Q]
        ######################
        z1=0
        z2=0
        for j in range(len(lis_Q)):
            F=[FET(lis_Q[j],self.b,z1,z2,Y[i]) for i in range(len(Y))]
            lis_F.append(F)
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.xlim(0,1.5*max(Y))
        plt.ylim(0,1.1*max(Y))
        yyy=posy*(dimy)
        xxx=posx*(1.5*max(Y))
        fac1=(dimy)/(self.y)
        yy=fac*self.y*fac1
        bb=fac*self.b*fac1
        array=np.array([[-bb/2+xxx,0+yyy],[-bb/2-z1*yy+xxx,yy+yyy],
                        [bb/2+z2*yy+xxx,yy+yyy],[bb/2+xxx,0+yyy]])
        shape=Polygon(array,color=lis_col[2])
        ax1.add_patch(shape)
        plt.plot([-bb/2-1.25*yy*z1+xxx,-bb/2+xxx,bb/2+xxx,bb/2+1.25*yy*z2+xxx],
                 [1.25*yy+yyy,0+yyy,0+yyy,1.25*yy+yyy],"-",color=lis_col[3])                    
        plt.plot(asintota,Y,"-",color="silver",lw=2)
    #liscol=["blue","black","magenta","red","orange","cyan","yellow"]
        for i in range(len(lis_Q)):
            FF=lis_F[i]
            plt.errorbar(FF,Y,linestyle='-',
                         color=lis_col[0],
                         label="F",lw=2)
        for i in range(len(lis_py)):
            lista=[lis_py[i],FET(self.Q,self.b,z1,z2,lis_py[i])]
            plt.plot([0,lista[1]],[lista[0],lista[0]],"--",color=lis_col[1])
            plt.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_col[1])
            plt.plot(lista[1],lista[0], "mo",markersize=11,mec="k",alpha=1,mew=1.5)
        plt. xlabel ("Fuerza específica (F)",style="italic",
                     fontweight="bold")
        plt. ylabel ("Tirante (Y)",style="italic",
                     fontweight="bold")
        plt.title("F vs Y",style="italic",
                  fontweight="bold")
        plt.legend(loc="upper right")
        #nx = 1
        plt.grid()
        plt.show()
    def show_EFY(self,alfa=1,posx=0.25,posy=0.78,tam=0.14,dimy=3*0.6):
    ###################################################
    ##########################
        y1=self.y
        #y2=self.results(decimal_number=6,tol=0.0001,show_iteration="FALSE")
        A=self.b*self.y
        V=self.Q/A
        g=9.807
        T=self.b
        F=(V)/(g*(A/T))**0.5
        y2=(self.y/2)*((1+8*F**2)**0.5-1)
        z1=0
        z2=0
        E1=EET(self.Q,self.b,z1,z2,alfa,y1)
        E2=EET(self.Q,self.b,z1,z2,alfa,y2)
        HT=E1-E2
        M1=FET(self.Q,self.b,z1,z2,y1)
        M2=FET(self.Q,self.b,z1,z2,y2)
###############################
        lis_pM=[[y1,M1],[y2,M2]]
        lis_pE=[[y1,EET(self.Q,self.b,z1,z2,alfa,y1)],
                [y2,EET(self.Q,self.b,z1,z2,alfa,y2)]]
        lis_colE=["orange","red","cyan","black"]
        lis_colM=["blue","red","cyan","black"]
        lis_datsec=[y1,self.b,z1,z2]
        lis_pos=[posx,posy,tam]
        #dimy=3*y1
        N=50
################################################
        inc=(dimy)/N
        Y=[inc*(i+1) for i in range(N)]
        asintota=[]
        for i in range(N):
            A=0.5*(z1+z2)*Y[i]**2+self.b*Y[i]
            ym=0.5*Y[i]
            asintota.append(ym*A)    
        F=[FET(self.Q,self.b,z1,z2,Y[i]) for i in range(len(Y))]
        E=[EET(self.Q,self.b,z1,z2,alfa,Y[i]) for i in range(len(Y))] 
##############################################
##########################################
        fig=plt.figure(constrained_layout=True,figsize=(10,5))    
        gs=GridSpec(10,18,figure=fig)
        fig1=fig.add_subplot(gs[1:9,1:9])
        fig2=fig.add_subplot(gs[1:9,9:17])
###########################################
        fig1.set_xlim(0,1.5*max(Y))
        fig1.set_ylim(0,1.1*max(Y))
        fig2.set_xlim(0,1.5*max(Y))
        fig2.set_ylim(0,1.1*max(Y))

        fig2.plot(asintota,Y, "-",color="silver",lw=2)
        fig2.plot(F,Y,linestyle ='-',color=lis_colM[0],lw=3,label="M")
        for i in range(len(lis_pM)):
            lista=lis_pM[i]
            fig2.plot([0,lista[1]],[lista[0],lista[0]],"--",color=lis_colM[1])
            fig2.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_colM[1])
            fig2.plot(lista[1],lista[0], "mo",markersize=11,mec="k",alpha=1,mew=1.5)
        fig2.set_title("F vs Y",style="italic",fontweight="bold")
        fig2.set_xlabel("Fuerza específica (F)",style="italic",fontweight="bold")
        fig2.grid(True)
########################################
#######energia específica
        y=lis_datsec[0]
        fac=lis_pos[2]
        yyy=lis_pos[1]*(dimy)
        xxx=lis_pos[0]*(1.5*max(Y))
        fac1=(dimy)/(y)
        y=fac*y*fac1
        b=fac*self.b*fac1
        array=np.array([[-b/2+xxx,0+yyy],[-b/2-z1*y+xxx,y+yyy],
                    [b/2+z2*y+xxx,y+yyy],[b/2+xxx,0+yyy]])
        shape=Polygon(array,color=lis_colE[2])
        fig1.add_patch(shape)
        fig1.plot([-b/2-1.25*y*z1+xxx,-b/2+xxx,b/2+xxx,b/2+1.25*y*z2+xxx],
                    [1.25*y+yyy,0+yyy,0+yyy,1.25*y+yyy],"-",color=lis_colE[3])                    
        fig1.plot(Y,Y,"-",color="silver",lw=2)
        fig1.plot(E,Y,"-",color=lis_colE[0],label="E",lw=3)
        for i in range(len(lis_pE)):
            lista=lis_pE[i]
            fig1.plot([lista[1],1.5*max(Y)],[lista[0],lista[0]],"--",color=lis_colE[1])
            fig1.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_colE[1])
            fig1.plot(lista[1],lista[0], "mo",markersize=11,mec="k",alpha=1,mew=1.5)
        fig1.plot([E1,E2],[y1-0.4,y1-0.4], "-",alpha=1,color="black",lw=1)
        fig1.plot([E2+0.5*(E1-E2),0.5*E1],[y1-0.4,y1-0.4+0.15], "-",alpha=1,color="black",lw=1)
        fig1.text(0.5*E1-0.3,y1-0.4+0.20,"hf[m]="+str(round(E1-E2,3)))
        #fig1.text(0.5*E1,y1-0.4+0.20,round(E1-E2,3))
        fig1.set_title("E vs Y",style="italic",fontweight="bold")
        fig1.set_ylabel("Tirante (Y)",style="italic",fontweight="bold")
        fig1.set_xlabel("Energía específica (E)",style="italic",fontweight="bold")
        fig2.grid(True)
        fig1.grid(True)
        fig1.legend(loc="upper right")
        fig2.legend(loc="upper right")
######################################################   
        plt.show()



##################################################
class Resalto_t:
    def __init__(self,Q,y,b,z1,z2):
        self.Q=Q
        self.y=y
        self.b=b
        self.z1=z1
        self.z2=z2
    def results(self,decimal_number=6,tol=0.0001,show_iteration="TRUE"):
        A=self.b*self.y+0.5*(self.z1+self.z2)*self.y**2
        V=self.Q/A
        g=9.81
        r=(V**2)/(2*g*self.y)
        z=(self.z1+self.z2)/2
        t=self.b/(z*self.y)
    ###########################
        error=500
        tol=0.0001
        J1=50
        n=10
        a=0
        b=500
        for i in range(5):
            xy=tabulador_ej(a,b,n,r,t)
            xx=xy[0]
            yy=xy[1]
            for j in range(len(xx)):
                if (yy[j]*yy[j+1]<0):
                    a=xx[j]
                    b=xx[j+1]
                    break
        J0=(a+b)/2 
        k=0
        if (show_iteration=="TRUE"):
            print("=======================Resalto hidráulico========================")
            print("N°".center(5),"J0".center(12),"E(J0)".center(12),"E'(J0)".center(12),
                  "J1".center(12),"error".center(12))
            print("=================================================================")
        while (error>tol):
            k=k+1
            EJ=E(J0,r,t)
            dEJ=dE(J0,r,t)
            J1=J0-EJ/dEJ
            error=abs(J1-J0)
            print(str(k+1).center(5),str(round(J0,decimal_number)).center(12),
                    str(round(EJ,decimal_number)).center(12),str(round(dEJ,decimal_number)).center(12),
                    str(round(J1,decimal_number)).center(12),str(round(error,decimal_number)).center(12))    
            J0=J1
        if (show_iteration=="TRUE"):
            print("=================================================================")
            print("Valor de J                  (J) :",J1)
            if (J1>1):
                print("Tirante supercrítico    (y1)[m] :",self.y)
                print("Tirante subcrítico      (y2)[m] :",self.y*J1)
            if (J1<1):
                print("Tirante supercrítico    (y1)[m] :",self.y*J1)
                print("Tirante subcrítico      (y2)[m] :",self.y)
        return self.y*J1        
    def show_FY(self,lis_py=[],lis_col=["blue","red","orange","gray"],posx=0.80,posy=0.5,tam=0.1,dimy=3*0.6):
        n=50
        fac=tam
        inc=(dimy)/n
        Y=[inc*(i+1) for i in range(n)]
        asintota=[]
        for i in range(n):
            A=0.5*(self.z1+self.z2)*Y[i]**2+self.b*Y[i]
            ym=(1/3+1/6*((self.b*Y[i])/A))*Y[i]
            asintota.append(ym*A)    
        lis_F=[]
        ######################
        lis_Q=[self.Q]
        ######################
        for j in range(len(lis_Q)):
            F=[FET(lis_Q[j],self.b,self.z1,self.z2,Y[i]) for i in range(len(Y))]
            lis_F.append(F)
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.xlim(0,1.5*max(Y))
        plt.ylim(0,1.1*max(Y))
        yyy=posy*(dimy)
        xxx=posx*(1.5*max(Y))
        fac1=(dimy)/(self.y)
        y=fac*self.y*fac1
        b=fac*self.b*fac1
        array=np.array([[-b/2+xxx,0+yyy],[-b/2-self.z1*y+xxx,y+yyy],
                        [b/2+self.z2*y+xxx,y+yyy],[b/2+xxx,0+yyy]])
        shape=Polygon(array,color=lis_col[2])
        ax1.add_patch(shape)
        plt.plot([-b/2-1.25*y*self.z1+xxx,-b/2+xxx,b/2+xxx,b/2+1.25*y*self.z2+xxx],
                 [1.25*y+yyy,0+yyy,0+yyy,1.25*y+yyy],"-",color=lis_col[3])                    
        plt.plot(asintota,Y,"-",color="silver",lw=2)
    #liscol=["blue","black","magenta","red","orange","cyan","yellow"]
        for i in range(len(lis_Q)):
            FF=lis_F[i]
            plt.errorbar(FF,Y,linestyle='-',
                         color=lis_col[0],
                         label="F",lw=2)
        for i in range(len(lis_py)):
            lista=[lis_py[i],FET(self.Q,self.b,self.z1,self.z2,lis_py[i])]
            plt.plot([0,lista[1]],[lista[0],lista[0]],"--",color=lis_col[1])
            plt.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_col[1])
            plt.plot(lista[1],lista[0], "mo",markersize=11,mec="k",alpha=1,mew=1.5)
        plt. xlabel ("Fuerza específica (F)",style="italic",
                     fontweight="bold")
        plt. ylabel ("Tirante (Y)",style="italic",
                     fontweight="bold")
        plt.title("F vs Y",style="italic",
                  fontweight="bold")
        plt.legend(loc="upper right")
        #nx = 1
        plt.grid()
        plt.show()
    def show_EFY(self,alfa=1,posx=0.25,posy=0.78,tam=0.14,dimy=3*0.6):
    ###################################################
    ##########################
        y1=self.y
        y2=self.results(decimal_number=6,tol=0.0001,show_iteration="FALSE")
        E1=EET(self.Q,self.b,self.z1,self.z2,alfa,y1)
        E2=EET(self.Q,self.b,self.z1,self.z2,alfa,y2)
        HT=E1-E2
        M1=FET(self.Q,self.b,self.z1,self.z2,y1)
        M2=FET(self.Q,self.b,self.z1,self.z2,y2)
###############################
        lis_pM=[[y1,M1],[y2,M2]]
        lis_pE=[[y1,EET(self.Q,self.b,self.z1,self.z2,alfa,y1)],[y2,EET(self.Q,self.b,self.z1,self.z2,alfa,y2)]]
        lis_colE=["orange","red","cyan","black"]
        lis_colM=["blue","red","cyan","black"]
        lis_datsec=[y1,self.b,self.z1,self.z2]
        lis_pos=[posx,posy,tam]
        #dimy=3*y1
        N=50
################################################
        inc=(dimy)/N
        Y=[inc*(i+1) for i in range(N)]
        asintota=[]
        for i in range(N):
            A=0.5*(self.z1+self.z2)*Y[i]**2+self.b*Y[i]
            ym=(1/3+1/6*((self.b*Y[i])/A))*Y[i]
            asintota.append(ym*A)    
        F=[FET(self.Q,self.b,self.z1,self.z2,Y[i]) for i in range(len(Y))]
        E=[EET(self.Q,self.b,self.z1,self.z2,alfa,Y[i]) for i in range(len(Y))] 
##############################################
##########################################
        fig=plt.figure(constrained_layout=True,figsize=(10,5))    
        gs=GridSpec(10,18,figure=fig)
        fig1=fig.add_subplot(gs[1:9,1:9])
        fig2=fig.add_subplot(gs[1:9,9:17])
###########################################
        fig1.set_xlim(0,1.5*max(Y))
        fig1.set_ylim(0,1.1*max(Y))
        fig2.set_xlim(0,1.5*max(Y))
        fig2.set_ylim(0,1.1*max(Y))

        fig2.plot(asintota,Y, "-",color="silver",lw=2)
        fig2.plot(F,Y,linestyle ='-',color=lis_colM[0],lw=3,label="M")
        for i in range(len(lis_pM)):
            lista=lis_pM[i]
            fig2.plot([0,lista[1]],[lista[0],lista[0]],"--",color=lis_colM[1])
            fig2.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_colM[1])
            fig2.plot(lista[1],lista[0], "mo",markersize=11,mec="k",alpha=1,mew=1.5)
        fig2.set_title("F vs Y",style="italic",fontweight="bold")
        fig2.set_xlabel("Fuerza específica (F)",style="italic",fontweight="bold")
        fig2.grid(True)
########################################
#######energia específica
        y=lis_datsec[0]
        fac=lis_pos[2]
        yyy=lis_pos[1]*(dimy)
        xxx=lis_pos[0]*(1.5*max(Y))
        fac1=(dimy)/(y)
        y=fac*y*fac1
        b=fac*self.b*fac1
        array=np.array([[-b/2+xxx,0+yyy],[-b/2-self.z1*y+xxx,y+yyy],
                    [b/2+self.z2*y+xxx,y+yyy],[b/2+xxx,0+yyy]])
        shape=Polygon(array,color=lis_colE[2])
        fig1.add_patch(shape)
        fig1.plot([-b/2-1.25*y*self.z1+xxx,-b/2+xxx,b/2+xxx,b/2+1.25*y*self.z2+xxx],
                    [1.25*y+yyy,0+yyy,0+yyy,1.25*y+yyy],"-",color=lis_colE[3])                    
        fig1.plot(Y,Y,"-",color="silver",lw=2)
        fig1.plot(E,Y,"-",color=lis_colE[0],label="E",lw=3)
        for i in range(len(lis_pE)):
            lista=lis_pE[i]
            fig1.plot([lista[1],1.5*max(Y)],[lista[0],lista[0]],"--",color=lis_colE[1])
            fig1.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_colE[1])
            fig1.plot(lista[1],lista[0], "mo",markersize=11,mec="k",alpha=1,mew=1.5)
        fig1.plot([E1,E2],[y1-0.4,y1-0.4], "-",alpha=1,color="black",lw=1)
        fig1.plot([E2+0.5*(E1-E2),0.5*E1],[y1-0.4,y1-0.4+0.15], "-",alpha=1,color="black",lw=1)
        fig1.text(0.5*E1-0.3,y1-0.4+0.20,"hf[m]="+str(round(E1-E2,3)))
        #fig1.text(0.5*E1,y1-0.4+0.20,round(E1-E2,3))
        fig1.set_title("E vs Y",style="italic",fontweight="bold")
        fig1.set_ylabel("Tirante (Y)",style="italic",fontweight="bold")
        fig1.set_xlabel("Energía específica (E)",style="italic",fontweight="bold")
        fig2.grid(True)
        fig1.grid(True)
        fig1.legend(loc="center right")
        fig2.legend(loc="center right")
######################################################   
        plt.show()




#Q=4#Caudal (Q)[m3/s]
#b=1#Ancho de solera (b)[m]
#z1=2#Talud 1
#z2=2#Talud 2
#y=0.3
#rr=Resalto_t(Q,y,b,z1,z2)
#y2=rr.results()
#rr.show_FY(lis_py=[y,y2],dimy=10*y)
#rr.show_EFY(dimy=5*y)
def N(y,d):
    A=0.25*d**2*acos(1-2*y/d)-0.5*d**2*(y/d-y**2/d**2)**0.5*(1-2*y/d)
    N=A/d**2
    return N
def K(y,d):
    K=1-0.5*(1/(y/d))+(2*(y/d)**0.5*(1-(y/d))**(3/2))/(3*N(y,d))
    return K
def E_RC(y,ydato,d,Q):
    g=9.807
    N1=N(ydato,d)
    N2=N(y,d)
    K1=K(ydato,d)
    K2=K(y,d)
    Ey=(K2*N1*N2*(y/ydato)-K1*N1**2)/((ydato/d)**4*(1-(N1)/(N2)))-(Q**2)/(g*ydato**5)
    return Ey 
def tabulador_erc(a,b,n,ydato,d,Q):
    inc=(b-a)/n
    x=[a+inc*i for i in range(n+1)]
    y=[E_RC(x[i],ydato,d,Q) for i in range(n+1)]
    return [x,y]

class Resalto_c:
    def __init__(self,Q,y,d):
        self.Q=Q
        self.y=y
        self.d=d
    def results(self,decimal_number=6,tol=0.0001,show_iteration="TRUE"):
        error=500
        n=5
        a=0.001
        b=self.d
        for i in range(5):
            xy=tabulador_erc(a,b,n,self.y,self.d,self.Q)
            xx=xy[0]
            yy=xy[1]
            for j in range(len(xx)):
                if (yy[j]*yy[j+1]<0):
                    a=xx[j]
                    b=xx[j+1]
                    break
        i=0
        if (show_iteration=="TRUE"):
            print("=====================Resalto hidráulico sección circular========================")
            print("N°".center(5),"a".center(12),"b".center(12),"ym".center(12),
                  "f(a)".center(12),"f(ym)".center(12),"error".center(12))
            print("================================================================================")
        lis_a=[]
        lis_b=[]
        while (error>tol):
            lis_a.append(a)
            lis_b.append(b)
            ym=0.5*(a+b)
            fa=E_RC(a,self.y,self.d,self.Q)
            fm=E_RC(ym,self.y,self.d,self.Q)
            if (fa*fm)<0:
                a=a
                b=ym
            else:
                a=ym
                b=b
            error=abs(ym-lis_a[i])
            if (show_iteration=="TRUE"):
                print(str(i).center(5),str(round(lis_a[i],decimal_number)).center(12),str(round(lis_b[i],decimal_number)).center(12),
                      str(round(ym,decimal_number)).center(12),str(round(fa,decimal_number)).center(12),
                      str(round(fm,decimal_number)).center(12),str(round(error,decimal_number)).center(12))    
            i=i+1
        if (show_iteration=="TRUE"):
            print("================================================================================")
        if (show_iteration=="TRUE"):
            if (self.y>ym):            
                print("Tirante supercrítico    (y1)[m] :",ym)
                print("Tirante subcrítico      (y2)[m] :",self.y)
            else:
                print("Tirante supercrítico    (y1)[m] :",self.y)
                print("Tirante subcrítico      (y2)[m] :",ym)
        return ym


#Q=1.5
#d=1.2
#y=0.6
#rr=Resalto_c(Q,y,d)
#rr.results()
















