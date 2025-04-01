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


class Intgraf:
    def __init__(self,Q,b,z1,z2,S0,n,y1,y2,nt):
        self.Q=Q
        self.b=b
        self.z1=z1
        self.z2=z2
        self.S0=S0
        self.n=n
        self.y1=y1
        self.y2=y2
        self.nt=nt
    def results(self,decimal_number=6,show_results="TRUE"):
        g=9.81
        deltay=abs(self.y1-self.y2)/self.nt
        if (self.y1<self.y2):
            Y=[self.y1+deltay*(i) for i in range(self.nt+1)]
        elif (self.y2<self.y1):
            Y=[self.y1-deltay*(i) for i in range(self.nt+1)]
        if (show_results=="TRUE"):
            print("y".center(12),"A[m]".center(12),"T[m]".center(12),
                  "R[m]".center(12),"v[m/s]".center(12),"Se".center(12),
                  "1-(Q^2T)/(gA^3)".center(12),"S0-Se".center(12),"f(y)".center(12),
                  "Δx=A".center(12),"x".center(12))
            print("===============================================================================================================================================")
        A=0.5*(self.z1+self.z2)*Y[0]**2+self.b*Y[0]
        P=((self.z1**2+1)**0.5+(self.z2**2+1)**0.5)*Y[0]+self.b
        R=A/P
        V=self.Q/A
        T=Y[0]*(self.z1+self.z2)+self.b
        Se=(V**2*self.n**2)/(R**(4/3))
        num=1-(self.Q**2*T)/(g*A**3)
        den=self.S0-Se
        fy=num/den
        rr=decimal_number
        if (show_results=="TRUE"):
            print(str(round(Y[0],rr)).center(12),str(round(A,rr)).center(12),str(round(T,rr)).center(12),
                  str(round(R,rr)).center(12),str(round(V,rr)).center(12),str(round(Se,rr)).center(12),
                  str(round(num,rr)).center(12),str(round(den,rr)).center(12),str(round(fy,rr)).center(12),"---".center(12),"---".center(12))
        lis_fy=[fy]
        X=[0]
        x=0
        for i in range(1,self.nt+1):
            A=0.5*(self.z1+self.z2)*Y[i]**2+self.b*Y[i]
            P=((self.z1**2+1)**0.5+(self.z2**2+1)**0.5)*Y[i]+self.b
            R=A/P
            V=self.Q/A
            T=Y[i]*(self.z1+self.z2)+self.b
            Se=(V**2*self.n**2)/(R**(4/3))
            num=1-(self.Q**2*T)/(g*A**3)
            den=self.S0-Se
            fy=num/den
            lis_fy.append(fy)
            deltax=deltay*(lis_fy[i-1]+lis_fy[i])/2
            x=x+deltax
            X.append(x)
            if (show_results=="TRUE"):
                print(str(round(Y[i],rr)).center(12),str(round(A,rr)).center(12),str(round(T,rr)).center(12),
                    str(round(R,rr)).center(12),str(round(V,rr)).center(12),str(round(Se,rr)).center(12),
                    str(round(num,rr)).center(12),str(round(den,rr)).center(12),str(round(fy,rr)).center(12),
                    str(round(deltax,rr)).center(12),str(round(x,rr)).center(12))
        #SS0=[abs(X[i])*S0 for i in range(len(X))]
        SS0=[0 for i in range(len(X))]
        ymax=max(X)*self.S0
        SSr=[ymax-X[i]*self.S0 for i in range(len(X))]
        #SSr=rev(SSr)
        YP=[Y[i]+SS0[i] for i in range(len(X))]
        YPr=[Y[i]+SSr[i] for i in range(len(X))]
        nn=15
        lisb=[]
        for i in range(nn):
            lisb.append((self.b/(nn-1))*(i))
        res1=makemat(0,nn,len(YP))
        res2=makemat(0,nn,len(YP))
        for i in range(nn):
            for j in range(len(YP)):
                res1[i][j]=YP[j]
                res2[i][j]=0
        return [X,YP,SS0,res1,res2,lisb,SSr,YPr]

    def show_Intgraf(self):
        resul=self.results(decimal_number=6,show_results="FALSE")
        X=resul[0]
        YP=resul[1]
        SS0=resul[2]
        res=resul[3]
        res2=resul[4]
        lisb=resul[5]
        SSr=resul[6]
        YPr=resul[7]
        fig=plt.figure(constrained_layout=True,figsize=(10,5))
        gs=GridSpec(20,18,figure=fig)
        ax=fig.add_subplot(gs[1:9,1:17])
        ax2=fig.add_subplot(gs[11:19,9:17])
        ax.plot(X,YPr,"-",color="blue",lw=4,label="agua")
        ax.plot(X,SSr,"-",color="silver",lw=3,label="pendiente")
        ax2.plot(X,YP,"-",color="orange",lw=4,label="agua")
        ax2.plot(X,SS0,"-",color="silver",lw=4,label="fondo")
        ax.set_xlabel("Distancia (x)[m]",style="italic",fontweight="bold")
        ax.set_ylabel("Altura (y)[m]",style="italic",fontweight="bold")
        ax.set_title("Integración gráfica",style="italic",fontweight="bold")
        ax2.set_xlabel("Distancia (x)[m]",style="italic",fontweight="bold")
        ax2.set_title("",style="italic",fontweight="bold")
        nx = 1
        ax.legend()
        ax2.legend()
        ax.grid(True)
        ax2.grid(True)
        ###########################
        datos=np.array(res)
        datos2 = np.array(res2)
        xi=np.array(X)
        yi=np.array(lisb)
        X,Y=np.meshgrid(xi,yi)
        ax=fig.add_subplot(2,2,3,projection='3d')
        ax.plot_wireframe(X,Y,datos)
        surf=ax.plot_surface(X,Y,datos,cmap='autumn_r')
        surf=ax.plot_surface(X,Y,datos2,cmap='Spectral_r')
        ax.set(zlim3d=(0,2*max(YP)), zlabel='y')
        ax.set_xlabel('(x)[m]')
        ax.set_ylabel('(b)[m]')
        ax.set_zlabel('yn')
        plt.show()

#Q=5#Caudal (Q)[m3/s]
#b=2.5#Ancho de solera (b)[m]
#z1=1.5#Talud 1
#z2=1.5#Talud 2
#S0=0.0005#Pendientel del fondo (S0)
#n=0.025#Coeficiente de Manning
#y1=1.4#Tirante inicial (y1)[m]
#y2=2.3#Tirante final (y2)[m]
#nt=10#Numero de tramos
################################
#res=Intgraf(Q,b,z1,z2,S0,n,y1,y2,nt)
#res.results()
#res.show_Intgraf()


class Directotramos:
    def __init__(self,Q,b,z1,z2,S0,n,y1,y2,nt):
        self.Q=Q
        self.b=b
        self.z1=z1
        self.z2=z2
        self.S0=S0
        self.n=n
        self.y1=y1
        self.y2=y2
        self.nt=nt
    def results(self,decimal_number=6,show_results="TRUE"):
        g=9.81
        dtex=8
        rr=decimal_number
        alfa=1
        deltay=abs(self.y1-self.y2)/self.nt
        if (self.y1<self.y2):
            Y=[self.y1+deltay*(i) for i in range(self.nt+1)]
        elif (self.y2<self.y1):
            Y=[self.y1-deltay*(i) for i in range(self.nt+1)]
        if (show_results=="TRUE"):
            print("y".center(dtex),"A[m]".center(dtex),"P[m]".center(dtex),
                  "R[m]".center(dtex),"R^(2/3)[m]".center(dtex),"v[m/s]".center(dtex),
                  "V^2/(2g)".center(dtex),"E".center(dtex),"ΔE".center(dtex),
                  "SE".center(dtex),"SEm".center(dtex),"S0-SEm".center(dtex),
                  "Δx[m]".center(dtex),"L[m]".center(dtex))
            print("===============================================================================================================================================")
        A=0.5*(self.z1+self.z2)*Y[0]**2+self.b*Y[0]
        P=((self.z1**2+1)**0.5+(self.z2**2+1)**0.5)*Y[0]+self.b
        R=A/P
        R23=R**(2/3)
        V=self.Q/A
        V2g=(V**2)/(2*g)
        E=Y[0]+alfa*(V**2)/(2*g)
        SE=(V**2*self.n**2)/(R**(4/3))
        print(str(round(Y[0],rr)).center(dtex),str(round(A,rr)).center(dtex),str(round(P,rr)).center(dtex),
             str(round(R,rr)).center(dtex),str(round(R23,rr)).center(dtex),str(round(V,rr)).center(dtex),
                str(round(V2g,rr)).center(dtex),str(round(E,rr)).center(dtex),"---".center(dtex),
                str(round(SE,rr)).center(dtex),"---".center(dtex),"---".center(dtex),
                "---".center(dtex),"---".center(dtex))
        lis_SE=[SE]
        lis_E=[E]
        lis_x=[0]
        X=[0]
        x=0
        for i in range(1,self.nt+1):
            A=0.5*(self.z1+self.z2)*Y[i]**2+self.b*Y[i]
            P=((self.z1**2+1)**0.5+(self.z2**2+1)**0.5)*Y[i]+self.b
            R=A/P
            R23=R**(2/3)
            V=self.Q/A
            V2g=(V**2)/(2*g)
            E=Y[i]+alfa*(V**2)/(2*g)
            lis_E.append(E)
            deltaE=lis_E[i]-lis_E[i-1]
            SE=(V**2*self.n**2)/(R**(4/3))
            lis_SE.append(SE)
            SEm=(lis_SE[i-1]+lis_SE[i])*0.5
            S0SE=self.S0-SEm
            deltax=(deltaE)/(S0SE)
            x=x+deltax
            X.append(x)
            if (show_results=="TRUE"):
                print(str(round(Y[i],rr)).center(dtex),str(round(A,rr)).center(dtex),str(round(P,rr)).center(dtex),
                    str(round(R,rr)).center(dtex),str(round(R23,rr)).center(dtex),str(round(V,rr)).center(dtex),
                    str(round(V2g,rr)).center(dtex),str(round(E,rr)).center(dtex),str(round(deltaE,rr)).center(dtex),
                    str(round(SE,rr)).center(dtex),str(round(SEm,rr)).center(dtex),str(round(S0SE,rr)).center(dtex),
                    str(round(deltax,rr)).center(dtex),str(round(x,rr)).center(dtex))
    
    #SS0=[abs(X[i])*S0 for i in range(len(X))]
        SS0=[0 for i in range(len(X))]
    #SS0=[0 for i in range(len(X))]
        ymax=max(X)*self.S0
        SSr=[ymax-X[i]*self.S0 for i in range(len(X))]
        YP=[Y[i]+SS0[i] for i in range(len(X))]
        YPr=[Y[i]+SSr[i] for i in range(len(X))]
        nn=15
        lisb=[]
        for i in range(nn):
            lisb.append((self.b/(nn-1))*(i))
        res1=makemat(0,nn,len(YP))
        res2=makemat(0,nn,len(YP))
        for i in range(nn):
            for j in range(len(YP)):
                res1[i][j]=YP[j]
                res2[i][j]=0       
    #return [X,YP,SS0,res1,res2,lisb]
        return [X,YP,SS0,res1,res2,lisb,SSr,YPr]

    def show_Directotramos(self):
        resul=self.results(decimal_number=6,show_results="FALSE")
        X=resul[0]
        YP=resul[1]
        SS0=resul[2]
        res=resul[3]
        res2=resul[4]
        lisb=resul[5]
        SSr=resul[6]
        YPr=resul[7]
        fig=plt.figure(constrained_layout=True,figsize=(10,5))
        gs=GridSpec(20,18,figure=fig)
        ax=fig.add_subplot(gs[1:9,1:17])
        ax2=fig.add_subplot(gs[11:19,9:17])
        ax.plot(X,YPr,"-",color="blue",lw=4,label="agua")
        ax.plot(X,SSr,"-",color="silver",lw=3,label="pendiente")
        ax2.plot(X,YP,"-",color="orange",lw=4,label="agua")
        ax2.plot(X,SS0,"-",color="silver",lw=4,label="fondo")
        ax.set_xlabel("Distancia (x)[m]",style="italic",fontweight="bold")
        ax.set_ylabel("Altura (y)[m]",style="italic",fontweight="bold")
        ax.set_title("Directo por tramos",style="italic",fontweight="bold")
        ax2.set_xlabel("Distancia (x)[m]",style="italic",fontweight="bold")
        ax2.set_title("",style="italic",fontweight="bold")
        nx = 1
        ax.legend()
        ax2.legend()
        ax.grid(True)
        ax2.grid(True)
        ###########################
        datos=np.array(res)
        datos2 = np.array(res2)
        xi=np.array(X)
        yi=np.array(lisb)
        X,Y=np.meshgrid(xi,yi)
        ax=fig.add_subplot(2,2,3,projection='3d')
        ax.plot_wireframe(X,Y,datos)
        surf=ax.plot_surface(X,Y,datos,cmap='autumn_r')
        surf=ax.plot_surface(X,Y,datos2,cmap='Spectral_r')
        ax.set(zlim3d=(0,2*max(YP)), zlabel='y')
        ax.set_xlabel('(x)[m]')
        ax.set_ylabel('(b)[m]')
        ax.set_zlabel('yn')
        plt.show()

#Q=1#Caudal (Q)[m3/s]
#b=0.8#Ancho de solera (b)[m]
#z1=1#Talud 1
#z2=1#Talud 2
#S0=0.01#Pendientel del fondo (S0)
#n=0.015#Coeficiente de Manning
#y1=0.447#Tirante inicial (y1)[m]
#y2=0.356#Tirante final (y2)[m]
#nt=10#Numero de tramos
#rr=Directotramos(Q,b,z1,z2,S0,n,y1,y2,nt)
#rr.results()
#rr.show_Directotramos()
def solve_tramosfijos(Q,S0,B,ZI,ZD,Y1,n,deltax):
    g=9.81
    A1=0.5*(ZI+ZD)*Y1**2+B*Y1
    P1=((ZI**2+1)**0.5+(ZD**2+1)**0.5)*Y1+B
    C=S0*deltax+Y1+(Q**2)/(2*g*A1**2)-((deltax*Q**2*n**2)/(2))*((P1**2)/(A1**5))**(2/3)
    error=500
    tol=0.0001
    Y2=0.59
    while (error>tol):
        A2=0.5*(ZI+ZD)*Y2**2+B*Y2
        P2=((ZI**2+1)**0.5+(ZD**2+1)**0.5)*Y2+B
        DA2=(ZI+ZD)*Y2+B
        DP2=(ZI**2+1)**0.5+(ZD**2+1)**0.5
        E=Y2+(Q**2)/(2*g*A2**2)+((deltax*Q**2*n**2)/(2))*((P2**2)/(A2**5))**(2/3)-C
        DE=1-(Q**2)/(g*A2**3)+((deltax*Q**2*n**2)/(2))*((4/3*A2**(10/3)*P2**(1/3)*DP2-10/3*P2**(4/3)*A2**(7/3)*DA2 )/(A2**(20/3)))
        Y3=Y2-E/DE
        error=abs(Y3-Y2)
        Y2=Y3
    return Y2
class Tramosfijos:
    def __init__(self,Q,S0,b,z1,z2,y1,n,deltax,nt):
        self.Q=Q
        self.S0=S0
        self.b=b
        self.z1=z1
        self.z2=z2
        self.y1=y1
        self.n=n
        self.deltax=deltax
        self.nt=nt
    def results(self,decimal_number=6,show_results="TRUE"):
        g=9.81
        dtex=12
        rr=decimal_number
        alfa=1
        X=[self.deltax*(i) for i in range(self.nt+1)]
        Y=[self.y1]
        for i in range(self.nt):
            Y.append(solve_tramosfijos(self.Q,self.S0,self.b,self.z1,self.z2,Y[i],self.n,self.deltax))
    #SS0=[abs(X[i])*S0 for i in range(nt+1)]
    #YP=[Y[i]+SS0[i] for i in range(nt+1)]
    #####################################################
        if (show_results=="TRUE"):    
            print("x".center(12),"y".center(12),"A[m2]".center(12),"P[m]".center(12),"R[m]".center(12),
                  "V[m/s]".center(12),"Se".center(12))
            print("=========================================================================================")
        for i in range(len(X)):
            A=0.5*(self.z1+self.z2)*Y[i]**2+self.b*Y[i]
            P=((self.z1**2+1)**0.5+(self.z2**2+1)**0.5)*Y[i]+self.b
            R=A/P
            V=self.Q/A
            Se=(V**2*self.n**2)/(R**(4/3))
            if (show_results=="TRUE"):
                print(str(X[i]).center(12),str(round(Y[i],rr)).center(12),
                      str(round(A,rr)).center(12),str(round(P,rr)).center(12),
                      str(round(R,rr)).center(12),str(round(V,rr)).center(12),
                      str(round(Se,rr)).center(12))        
        
    #SS0=[abs(X[i])*S0 for i in range(len(X))]
        SS0=[0 for i in range(len(X))]
    #YP=[Y[i]+SS0[i] for i in range(len(X))]
        ymax=max(X)*self.S0
        SSr=[ymax-X[i]*self.S0 for i in range(len(X))]
        YP=[Y[i]+SS0[i] for i in range(len(X))]
        YPr=[Y[i]+SSr[i] for i in range(len(X))]    
        nn=15
        lisb=[]
        for i in range(nn):
            lisb.append((self.b/(nn-1))*(i))
        res1=makemat(0,nn,len(YP))
        res2=makemat(0,nn,len(YP))
        for i in range(nn):
            for j in range(len(YP)):
                res1[i][j]=YP[j]
                res2[i][j]=0       
        return [X,YP,SS0,res1,res2,lisb,SSr,YPr]
    def show_Tramosfijos(self):
        resul=self.results(decimal_number=6,show_results="FALSE")
        X=resul[0]
        YP=resul[1]
        SS0=resul[2]
        res=resul[3]
        res2=resul[4]
        lisb=resul[5]
        SSr=resul[6]
        YPr=resul[7]
        fig=plt.figure(constrained_layout=True,figsize=(10,5))
        gs=GridSpec(20,18,figure=fig)
        ax=fig.add_subplot(gs[1:9,1:17])
        ax2=fig.add_subplot(gs[11:19,9:17])
        ax.plot(X,YPr,"-",color="blue",lw=4,label="agua")
        ax.plot(X,SSr,"-",color="silver",lw=3,label="pendiente")
        ax2.plot(X,YP,"-",color="orange",lw=4,label="agua")
        ax2.plot(X,SS0,"-",color="silver",lw=4,label="fondo")
        ax.set_xlabel("Distancia (x)[m]",style="italic",fontweight="bold")
        ax.set_ylabel("Altura (y)[m]",style="italic",fontweight="bold")
        ax.set_title("Tramos fijos",style="italic",fontweight="bold")
        ax2.set_xlabel("Distancia (x)[m]",style="italic",fontweight="bold")
        ax2.set_title("",style="italic",fontweight="bold")
        nx = 1
        ax.legend()
        ax2.legend()
        ax.grid(True)
        ax2.grid(True)
        ###########################
        datos=np.array(res)
        datos2 = np.array(res2)
        xi=np.array(X)
        yi=np.array(lisb)
        X,Y=np.meshgrid(xi,yi)
        ax=fig.add_subplot(2,2,3,projection='3d')
        ax.plot_wireframe(X,Y,datos)
        surf=ax.plot_surface(X,Y,datos,cmap='autumn_r')
        surf=ax.plot_surface(X,Y,datos2,cmap='Spectral_r')
        ax.set(zlim3d=(0,2*max(YP)), zlabel='y')
        ax.set_xlabel('(x)[m]')
        ax.set_ylabel('(b)[m]')
        ax.set_zlabel('yn')
        plt.show()

#Q=2#Caudal (Q)[m3/s]
#b=1#Ancho de solera (b)[m]
#z1=2#Talud izquierdo
#z2=2#Talud derecho
#S0=0.0005#Pendientel del fondo (S0)
#n=0.025#Coeficiente de Manning
#y1=1.5#Tirante inicial (y1)[m]
#nt=10#Numero de tramos
#deltax=-200 #Incremento de x
################################################
#rr=Tramosfijos(Q,S0,b,z1,z2,y1,n,deltax,nt)
#rr.results()
#rr.show_Tramosfijos()

def phiz(z):
    zz=[-50,-30,-20,-15,-12,-10,-9,-8,-7,-6.5,-6,-5.5,-5,-4.8,-4.6,-4.4,-4.2,
        -4,-3.8,-3.6,-3.4,-3.2,-3,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,
        -2,-1.95,-1.90,-1.85,-1.80,-1.75,-1.70,-1.65,-1.60,-1.55,-1.50,
        -1.48,-1.46,-1.44,-1.42,-1.40,
        -1.38,-1.36,-1.34,-1.32,-1.30,-1.28,-1.26,-1.24,-1.22,-1.20,-1.18,-1.16,
        -1.14,-1.12,-1.10,-1.08,-1.06,-1.04,-1.02,-1.00,-0.95,-0.90,
        -0.85,-0.80,-0.75,-0.70,-0.65,-0.60,-0.55,-0.50,-0.45,-0.40,-0.35,
        -0.30,-0.25,-0.20,-0.15,-0.10,-0,
        0,0.1,0.2,0.25,0.30,0.35,0.40,0.45,0.50,0.52,0.54,0.56,0.58,0.60,0.62,
        0.64,0.65,0.68,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,
        0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.905,0.910,0.915,
        0.920,0.925,0.930,0.935,0.940,0.945,0.950,0.952,0.954,0.956,0.958,0.960,
        0.962,0.964,0.966,0.968,0.970,0.971,0.972,0.973,0.974,0.975,0.976,
        0.977,0.978,0.979,0.980,0.981,0.982,0.983,0.984,0.985,0.986,0.987,
        0.988,0.989,0.990,0.991,0.992,0.993,0.994,0.995,0.996,0.997,0.998,
        0.999,
        1.001,1.002,1.003,1.004,1.005,1.006,1.007,1.008,1.009,1.010,1.011,
        1.012,1.013,1.014,1.015,1.016,1.017,1.018,1.019,1.020,1.021,1.022,
        1.023,1.024,1.025,1.026,1.027,1.028,1.029,1.030,1.031,1.032,1.033,
        1.034,1.035,1.036,1.037,1.038,1.039,1.040,1.041,1.042,1.043,1.044,
        1.045,1.046,1.047,1.048,1.049,1.050,1.052,1.054,1.056,1.058,1.060,
        1.062,1.064,1.066,1.068,1.070,1.072,1.074,1.076,1.078,1.080,1.082,
        1.084,1.086,1.088,1.090,1.092,1.094,1.096,1.098,1.1,1.105,1.110,
        1.115,1.120,1.125,1.130,1.135,1.140,1.145,1.150,1.155,1.160,1.165,
        1.170,1.175,1.180,1.185,1.190,1.195,1.200,1.21,1.22,1.23,1.24,1.25,
        1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,
        1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.52,
        1.54,1.56,1.58,1.60,1.62,1.64,1.66,1.68,1.70,1.72,1.74,1.76,1.78,
        1.80,1.82,1.84,1.86,1.88,1.90,1.92,1.94,1.96,1.98,2,2.05,2.10,2.15,
        2.20,2.25,2.30,2.35,2.40,2.45,2.50,2.55,2.60,2.65,2.70,2.75,2.80,
        2.85,2.90,2.95,3.00,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,
        4.3,4.4,4.5,4.6,4.7,4.8,4.9,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,12,15,
        20.2,30,50,100]
    phiz=[0.0002,0.0006,0.0013,0.0022,0.0035,0.0050,0.0062,0.0078,0.0102,
          0.0118,0.0139,0.0165,0.0199,0.0216,0.0235,0.0257,0.0282,0.0311,
          0.0344,0.0383,0.0428,0.0482,0.0548,0.0585,0.0626,0.0672,
          0.0723,0.0780,0.0845,0.0916,0.0996,0.1088,0.1192,0.1249,0.1311,
          0.1377,0.1147,0.1523,0.1605,0.1692,0.1787,0.1889,0.1999,0.2045,
          0.2093,0.2143,0.2194,0.2246,0.2301,0.2357,0.2414,0.2474,0.2536,
          0.2599,0.2665,0.2733,0.2802,0.2875,0.2949,0.3026,0.3105,0.3187,
          0.3272,0.3359,0.3449,0.3541,0.3637,0.3736,0.3995,0.4274,0.4574,
          0.4894,0.5234,0.5597,0.5979,0.6381,0.6801,0.7238,0.7689,0.8154,
          0.8629,0.9112,0.9603,1.0096,1.0593,1.1092,1.2092,
          0,0.1,0.2004,0.2510,0.3021,0.3538,0.4066,0.4608,0.5168,0.5399,
          0.5634,0.5874,0.6120,0.6371,0.6630,0.6897,0.7173,0.7459,0.7757,
          0.7910,0.8068,0.8230,0.8396,0.8566,0.8742,0.8923,0.9110,
          0.9304,0.9505,0.9714,0.9932,1.0160,1.0399,1.0651,1.0918,
          1.1202,1.1505,1.1831,1.2184,1.2373,1.2571,1.2779,1.2999,1.3232,1.3479,
          1.3744,1.4025,1.4336,1.4670,1.4813,1.4962,1.5117,1.5279,1.5448,
          1.5626,1.5813,1.6011,1.6220,1.6442,1.6558,1.6678,1.6803,1.6932,
          1.7066,1.7206,1.7351,1.7503,1.7661,1.7827,1.8001,1.8185,1.8379,1.8584,
          1.8803,1.9036,1.9287,1.9557,1.9850,2.0171,2.0526,2.0922,2.1370,
          2.1887,2.2498,2.3246,2.4208,2.5563,2.7877,
          2.1837,1.9530,1.8182,1.7226,1.6486,1.5881,1.5371,1.4929,1.4540,
          1.4192,1.3878,1.3591,1.3327,1.3083,1.2857,1.2645,1.2446,1.2259,1.2082,
          1.1914,1.1755,1.1603,1.1458,1.1320,1.1187,1.1060,1.0937,
          1.0819,1.0706,1.0596,1.0490,1.0387,1.0288,1.0191,1.0098,1.0007,
          0.9919,0.9634,0.9750,0.9669,0.9590,0.9513,0.9438,0.9354,0.9293,
          0.9223,0.9154,0.9087,0.9022,0.8958,0.8834,0.8714,0.8599,0.8499,0.8382,
          0.8279,0.8180,0.8084,0.7990,0.7900,0.7813,0.7728,0.7645,0.7565,
          0.7487,0.7411,0.7337,0.7265,0.7194,0.7126,0.7059,0.6993,0.6929,
          0.6867,0.6806,0.6659,0.6519,0.6387,0.6260,0.6139,0.6025,0.5913,
          0.5808,0.5707,0.5608,0.5514,0.5423,0.5335,0.5251,0.5169,0.5090,
          0.5014,0.4939,0.4868,0.4798,0.4664,0.4538,0.4419,0.4306,0.4196,
          0.4096,0.3998,0.3905,0.3816,0.3731,0.3649,0.3570,0.3495,
          0.3422,0.3352,0.3285,0.3220,0.3158,0.3098,0.3039,0.2983,
          0.2928,0.2875,0.2824,0.2775,0.2680,0.2727,0.2635,0.2591,0.2548,
          0.2466,0.2389,0.2315,0.2246,0.2179,0.2116,0.2056,0.1999,0.1944,
          0.1892,0.1842,0.1794,0.1748,0.1704,0.1662,0.1621,0.1582,0.1545,
          0.1509,0.1474,0.1440,0.1408,0.1377,0.1347,0.1318,0.1249,0.1186,
          0.1128,0.1074,0.1024,0.0978,0.0935,0.0894,0.0857,0.0821,0.0788,
          0.0757,0.0728,0.0700,0.0674,0.0650,0.0626,0.0604,0.0584,0.0564,0.527,
          0.0494,0.0464,0.0437,0.0412,0.0389,0.0368,0.0349,0.0331,0.0315,
          0.0299,0.0285,0.0272,0.0259,0.0248,0.0237,0.0227,0.0218,0.0209,
          0.0201,0.0166,0.0139,0.0118,0.0102,0.0089,0.0077,0.0069,0.0062,
          0.0055,0.0050,0.0035,0.0022,0.0013,0.0006,0.0002,0.0001]
    res=interl(zz,phiz,z)
    return res
def interl(x,y,dx):
    for i in range(len(x)-1):
        if (x[i]<=dx<=x[i+1]):
            yy=((y[i+1]-y[i])/(x[i+1]-x[i]))*(dx-x[i])+y[i]
            break
    return yy
class Bresse:
    def __init__(self,Q,b,S0,n,yn,y1,y2,nt):
        self.Q=Q
        self.b=b
        self.S0=S0
        self.n=n
        self.yn=yn
        self.y1=y1
        self.y2=y2
        self.nt=nt
    def results(self,decimal_number=6,show_results="TRUE"):
        g=9.81
        dtex=12
        rr=decimal_number
        alfa=1
        deltay=abs(self.y1-self.y2)/self.nt
        if (self.y1<self.y2):
            Y=[self.y1+deltay*(i) for i in range(self.nt+1)]
        elif (y2<y1):
            Y=[self.y1-deltay*(i) for i in range(self.nt+1)]
        if (show_results=="TRUE"):
            print("y".center(dtex),"Z=y/yn".center(dtex),
                  "φ(Z)".center(dtex),"deltax".center(dtex),"x".center(dtex))
            print("===================================================================")
        y=0.5*(self.y1+self.y2)
        C=(y**(1/6))/(self.n)
        Z=Y[0]/self.yn
        FZ=phiz(Z)
        deltax=(self.yn)/(self.S0)*Z-self.yn*(1/self.S0-(C**2)/g)*FZ
        deltaxF=deltax
        print(str(round(Y[0],rr)).center(dtex),str(round(Z,rr)).center(dtex),
              str(round(FZ,rr)).center(dtex),str(round(deltax,rr)).center(dtex),"0".center(dtex)) 
        X=[0]
        x=0
        for i in range(1,self.nt+1):
            y=0.5*(self.y1+self.y2)
            C=(y**(1/6))/(self.n)
            Z=Y[i]/self.yn
            FZ=phiz(Z)
            deltax=(self.yn)/(self.S0)*Z-self.yn*(1/self.S0-(C**2)/g)*FZ
            x=abs(deltaxF-deltax)
            X.append(x)
            if (show_results=="TRUE"):
                print(str(round(Y[i],rr)).center(dtex),str(round(Z,rr)).center(dtex),
                    str(round(FZ,rr)).center(dtex),str(round(deltax,rr)).center(dtex),str(round(x,rr)).center(dtex))
                
    #SS0=[abs(X[i])*S0 for i in range(len(X))]
        SS0=[0 for i in range(len(X))]
    #YP=[Y[i]+SS0[i] for i in range(len(X))]
    ##########################################
        ymax=max(X)*self.S0
        SSr=[ymax-X[i]*self.S0 for i in range(len(X))]
        YP=[Y[i]+SS0[i] for i in range(len(X))]
        YPr=[Y[i]+SSr[i] for i in range(len(X))]
    ####################
        nn=15
        lisb=[]
        for i in range(nn):
            lisb.append((self.b/(nn-1))*(i))
        res1=makemat(0,nn,len(YP))
        res2=makemat(0,nn,len(YP))
        for i in range(nn):
            for j in range(len(YP)):
                res1[i][j]=YP[j]
                res2[i][j]=0       
        return [X,YP,SS0,res1,res2,lisb,SSr,YPr]
    def show_Bresse(self):
        resul=self.results(decimal_number=6,show_results="FALSE")
        X=resul[0]
        YP=resul[1]
        SS0=resul[2]
        res=resul[3]
        res2=resul[4]
        lisb=resul[5]
        SSr=resul[6]
        YPr=resul[7]
        fig=plt.figure(constrained_layout=True,figsize=(10,5))
        gs=GridSpec(20,18,figure=fig)
        ax=fig.add_subplot(gs[1:9,1:17])
        ax2=fig.add_subplot(gs[11:19,9:17])
        ax.plot(X,YPr,"-",color="blue",lw=4,label="agua")
        ax.plot(X,SSr,"-",color="silver",lw=3,label="pendiente")
        ax2.plot(X,YP,"-",color="orange",lw=4,label="agua")
        ax2.plot(X,SS0,"-",color="silver",lw=4,label="fondo")
        ax.set_xlabel("Distancia (x)[m]",style="italic",fontweight="bold")
        ax.set_ylabel("Altura (y)[m]",style="italic",fontweight="bold")
        ax.set_title("Bresse",style="italic",fontweight="bold")
        ax2.set_xlabel("Distancia (x)[m]",style="italic",fontweight="bold")
        ax2.set_title("",style="italic",fontweight="bold")
        nx = 1
        ax.legend()
        ax2.legend()
        ax.grid(True)
        ax2.grid(True)
        ###########################
        datos=np.array(res)
        datos2 = np.array(res2)
        xi=np.array(X)
        yi=np.array(lisb)
        X,Y=np.meshgrid(xi,yi)
        ax=fig.add_subplot(2,2,3,projection='3d')
        ax.plot_wireframe(X,Y,datos)
        surf=ax.plot_surface(X,Y,datos,cmap='autumn_r')
        surf=ax.plot_surface(X,Y,datos2,cmap='Spectral_r')
        ax.set(zlim3d=(0,2*max(YP)), zlabel='y')
        ax.set_xlabel('(x)[m]')
        ax.set_ylabel('(b)[m]')
        ax.set_zlabel('yn')
        plt.show()



#Q=10 #Caudal (Q)[m3/s]
#b=10 #Ancho de solera (b)[m]
#S0=0.0004 #Pendientel del fondo (S0)
#n=0.030 #Coeficiente de Manning
#yn=1.409 #Tirante normal (yn)[m]
#y1=1.420 #Tirante inicial (y1)[m]
#y2=3 #Tirante final (y2)[m]
#nt=10 #Numero de tramos
################################################
#rr=Bresse(Q,b,S0,n,yn,y1,y2,nt)
#rr.results()
#rr.show_Bresse()







