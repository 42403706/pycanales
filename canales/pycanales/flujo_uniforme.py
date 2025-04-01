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
def interl(x,y,dx):
    for i in range(len(x)-1):
        if (x[i]<=dx<=x[i+1]):
            yy=((y[i+1]-y[i])/(x[i+1]-x[i]))*(dx-x[i])+y[i]
            break
    return yy
def epsilons(z,by):
    lby=[0,1,1.5,2,4,6,8,10]
    epsz0=[0.2,0.46,0.59,0.67,0.73,0.73,0.73,0.73]
    epsz15=[0.56,0.68,0.70,0.710,0.715,0.720,0.730,0.740]
    epsz20=[0.64,0.71,0.72,0.73,0.74,0.75,0.76,0.77]
    if (z==0):
        eps=interl(lby,epsz0,by)
    elif (z==1.5):
        eps=interl(lby,epsz15,by)
    elif (z==2):
        eps=interl(lby,epsz20,by)    
    return eps


def epsilonp(z,by):
    lby=[0,1,1.5,2,3,4,6,8,10]
    epsz0=[0.1,0.41,0.53,0.65,0.84,0.9,0.9,0.9,0.9]
    epsz15=[0.3,0.79,0.84,0.9,0.92,0.96,0.98,0.99,1]
    if (z==0):
        epp=interl(lby,epsz0,by)
    elif (z==1.5 or z==2 ):
        epp=interl(lby,epsz15,by)    
    return epp
def newton(f,df,x0,tol,de):
    error=500
    i=0
    print("N0".center(5),"x0".center(12),"f(x0)".center(12),"f'(x0)".center(12),
          "x1".center(12),"error".center(12))
    print("=========================================================================")
    while (error>tol):
        i=i+1
        if (i>20):
            break
        x1=x0-(f(x0))/(df(x0))
        error=abs(x1-x0)
        print(str(i).center(5),str(round(x0,de)).center(12),
              str(round(f(x0),de)).center(12),
              str(round(df(x0),de)).center(12),
              str(round(x1,de)).center(12),str(round(error,de)).center(12))     
        x0=x1
    print("=========================================================================")    
    return x1
def caudal_sc(y,d,S,n):
    if (y==0):
        Q=0
    else:
        A=0.25*d**2*acos(1-2*y/d)-0.5*d**2*(y/d-y**2/d**2)**0.5*(1-2*y/d)
        P=d*acos(1-2*y/d)
        Q=(A**(5/3)*S**(1/2))/(P**(2/3)*n)
    return Q
def velocidad_sc(y,d,S,n):
    if (y==0):
        V=0
    else:
        A=0.25*d**2*acos(1-2*y/d)-0.5*d**2*(y/d-y**2/d**2)**0.5*(1-2*y/d)
        P=d*acos(1-2*y/d)
        V=(A**(2/3)*S**(1/2))/(P**(2/3)*n)
    return V
def lisQY(d,n,S,N):
    inc=(d)/N
    Y=[inc*(i) for i in range(N+1)]
    QQ=[0]
    for i in range(1,N+1):            
        QQ.append(caudal_sc(Y[i],d,S,n))
    return [QQ,Y]
def lisVY(d,n,S,N):
    inc=(d)/N
    Y=[inc*(i) for i in range(N+1)]
    VV=[0]
    for i in range(1,N+1):
        VV.append(velocidad_sc(Y[i],d,S,n))
    return [VV,Y]
#####################################    
class Tiranten_ty:
    def __init__(self,Q,b,z1,z2,n,sf):
        self.Q=Q
        self.b=b
        self.z1=z1
        self.z2=z2
        self.n=n
        self.sf=sf        
    def results(self,y0=0.5,decimal_number=6,tol=0.0001,show_iteration="FALSE",show_results="TRUE"):
        error=50
        if (show_iteration=="TRUE"):
            print("yi[m]".center(10),"A[m2]".center(10),
                  "dA".center(10),"P".center(10),"dP".center(10),
                  "f(y)".center(10),"df(y)".center(10),"yi+1[m]".center(10),
                  "error".center(10))
            print(" ==================================================================================================")
        lis_y0=[]
        lis_A=[]
        lis_dA=[]
        lis_P=[]
        lis_dP=[]
        lis_fy=[]
        lis_dfy=[]
        lis_y1=[]
        lis_error=[]
        while (error>tol):
            A=self.b*y0+0.5*(self.z1+self.z2)*y0**2
            P=((self.z1**2+1)**0.5+(self.z2**2+1)**0.5)*y0+self.b
            dA=self.b+(self.z1+self.z2)*y0
            dP=(self.z1**2+1)**0.5+(self.z2**2+1)**0.5
            fy=(A**(5/3))/(P**(2/3))-(self.Q*self.n)/(self.sf**(1/2))
            dfy=(5/3)*A**(2/3)*P**(-2/3)*dA-(2/3)*A**(5/3)*P**(-5/3)*dP
            y1=y0-fy/dfy
            error=abs(y1-y0)
            lis_y0.append(y0)
            lis_A.append(A)
            lis_dA.append(dA)
            lis_P.append(P)
            lis_dP.append(dP)
            lis_fy.append(fy)
            lis_dfy.append(dfy)
            lis_y1.append(y1)
            lis_error.append(error)
            if (show_iteration=="TRUE"):
                print(str(round(y0,decimal_number)).center(10),
                      str(round(A,decimal_number)).center(10),
                      str(round(dA,decimal_number)).center(10),str(round(P,decimal_number)).center(10),
                      str(round(dP,decimal_number)).center(10),
                      str(round(fy,decimal_number)).center(10),str(round(dfy,decimal_number)).center(10),
                      str(round(y1,decimal_number)).center(10),
                      str(round(error,decimal_number)).center(10))    
            y0=y1
        if (show_iteration=="TRUE"):
            print(" ==================================================================================================")       
        ff=len(lis_y0)
        mat_ite=makemat(0,ff,9)
        for i in range(ff):
            mat_ite[i][0]=lis_y0
            mat_ite[i][1]=lis_A
            mat_ite[i][2]=lis_dA
            mat_ite[i][3]=lis_P
            mat_ite[i][4]=lis_dP
            mat_ite[i][5]=lis_fy
            mat_ite[i][6]=lis_dfy
            mat_ite[i][7]=lis_y1
            mat_ite[i][8]=lis_error

        g=9.807
        T=self.z1*y0+self.b+self.z2*y0
        R=A/P
        V=self.Q/A
        F=(V)/((g*(A)/(T))**0.5)
        E=y0+(V**2)/(2*g)
        if (show_results=="TRUE"):
            print(" Tirante normal          yn[m] :",y0)
            print(" Área hidráulica         A[m2] :",A)
            print(" Espejo de agua           T[m] :",T)
            print(" Número de Froude            F :",F)
            print(" Perímetro mojado         P[m] :",P)
            print(" Radio hidráulico         R[m] :",R)
            print(" Velocidad              V[m/s] :",V)
            print(" Energía específica E[m-kg/kg] :",E)
        return [y1,A,T,F,P,R,V,E,mat_ite]
    def show_section3d(self,lis_col=["red","lightcoral","cyan","blue"]):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        L=10
        rr=0.25
        h=y+rr*y
        points = np.array([[-0.5*L,-0.5*self.b,0],
                          [0.5*L, -0.5*self.b, 0 ],
                          [0.5*L,-0.5*self.b-self.z1*h, h],
                          [-0.5*L,-0.5*self.b-self.z1*h,h],
                          [-0.5*L,0.5*self.b,0],
                          [0.5*L, 0.5*self.b, 0 ],
                          [0.5*L,0.5*self.b+self.z2*h, h],
                          [-0.5*L,0.5*self.b+self.z2*h,h] 
                           ])
        #ancho=(z1*h+b+z2*h)/15
        ancho=0.2*(self.z1*y+self.b+self.z2*y)
        anchov=0.2*(y)
        points2 = np.array([[-0.5*L,-0.5*self.b-ancho-self.z1*ancho,-anchov],
                          [0.5*L, -0.5*self.b-ancho-self.z1*ancho, -anchov],
                          [0.5*L,-0.5*self.b-self.z1*h-ancho, h],
                          [-0.5*L,-0.5*self.b-self.z1*h-ancho,h],
                          [-0.5*L,0.5*self.b+ancho+self.z2*ancho,-anchov],
                          [0.5*L, 0.5*self.b+ancho+self.z2*ancho, -anchov],
                          [0.5*L,0.5*self.b+self.z2*h+ancho, h],
                          [-0.5*L,0.5*self.b+self.z2*h+ancho,h] 
                           ])
    
        points1 = np.array([[-0.5*L,-0.5*self.b,0],
                          [0.5*L, -0.5*self.b, 0 ],
                          [0.5*L,-0.5*self.b-self.z1*y, y],
                          [-0.5*L,-0.5*self.b-self.z1*y,y],
                          [-0.5*L,0.5*self.b,0],
                          [0.5*L, 0.5*self.b, 0 ],
                          [0.5*L,0.5*self.b+self.z2*y, y],
                          [-0.5*L,0.5*self.b+self.z2*y,y] 
                           ])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        r = [-5,5]
        X, Y = np.meshgrid(r,r)
        ax.set(xlim3d=(-0.5*L,0.5*L), xlabel='L')
        ax.set(ylim3d=(-0.5*self.b-self.z1*h-ancho-0.1*h,0.5*self.b+self.z2*h+ancho+0.1*h ),
               ylabel='')
        ax.set(zlim3d=(-ancho-0.3*h,h+0.3*h), zlabel='y')
        # vértices 
        #ax.scatter3D(points[:,0],points[:,1],points[:,2]) 
        #Lista de polígonos de los lados de la figura
        verts=[[points[0],points[1],points[2],points[3]],
               [points[4],points[5],points[6],points[7]],
               [points[0],points[1],points[5],points[4]]
               ]
        verts2=[[points2[0],points2[1],points2[2],points2[3]],
               [points2[4],points2[5],points2[6],points2[7]],
               [points2[0],points2[1],points2[5],points2[4]],
               [points[2],points[3],points2[3],points2[2]],
               [points[6],points[7],points2[7],points2[6]]  
               ]
        verts1=[[points1[0],points1[1],points1[2],points1[3]],
                [points1[4],points1[5],points1[6],points1[7]],
                [points1[0],points1[1],points1[5],points1[4]],
                [points1[0],points1[4],points1[7],points1[3]],
                [points1[1],points1[5],points1[6],points1[2]],
                [points1[2],points1[6],points1[7],points1[3]]
                ]
    #lados de la trama
        ax.add_collection3d(Poly3DCollection(verts,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))

        ax.add_collection3d(Poly3DCollection(verts2,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))
        ax.add_collection3d(Poly3DCollection(verts1,
                                             facecolors=lis_col[2],linewidths=1,
                                             edgecolors=lis_col[3], alpha=0.5))
        if (self.b==0):
            plt.title("Tirante normal sección triangular")
        elif (self.z1==0 and self.z2==0):
            plt.title("Tirante normal sección rectangular")
        elif (self.z1!=0 or self.z2!=0):
            plt.title("Tirante normal sección trapezoidal")    
        plt.show()
    def show_section(self,outline_color="silver",fill_color="cyan",width_size=7,height_size=5):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        T=self.z1*y+self.b+self.z2*y
        fig1=plt.figure(figsize=(width_size,height_size)) 
        ax1=fig1.add_subplot(xticks=[], yticks=[])
        dx1=-self.b/2-1.25*y*self.z1 
        dx2=+self.b/2+1.25*y*self.z2
        plt.xlim(dx1-0.25*T,dx2+0.25*T)
        plt.ylim(-0.7*T+0.5*y,0.7*T+0.5*y)
        yyy=0
        xxx=0
        array=np.array([[-self.b/2+xxx,0+yyy],[-self.b/2-self.z1*y+xxx,y+yyy],
                          [self.b/2+self.z2*y+xxx,y+yyy],[self.b/2+xxx,0+yyy]])
        shape=Polygon(array,color=fill_color)
        ax1.add_patch(shape)
        plt.plot([-self.b/2-1.25*y*self.z1+xxx,-self.b/2+xxx,self.b/2+xxx,self.b/2+1.25*y*self.z2+xxx],
                   [1.25*y+yyy,0+yyy,0+yyy,1.25*y+yyy],"-",color=outline_color,lw=6)
        dy=0.05*abs(-dx1+dx2)
        if (self.b!=0):
            plt.plot([-self.b/2,self.b/2],[-dy,-dy],"-",alpha=1,color="black",lw=1)
            plt.plot([-self.b/2,-self.b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.plot([self.b/2,self.b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.text(0,-3*dy,str(round(self.b,3))+"  m",style="italic")
            plt.text(0-0.25*self.b,-3*dy,"b=",style="italic")

        plt.plot([self.b/2+self.z2*y+dy,self.b/2+self.z2*y+dy],[0,y],"-",alpha=1,color="black",lw=1)
        plt.plot([self.b/2+self.z2*y+dy-0.3*dy ,self.b/2+self.z2*y+dy+0.3*dy],[y,y],"-",alpha=1,color="black",lw=1)
        plt.plot([self.b/2+self.z2*y+dy-0.3*dy ,self.b/2+self.z2*y+dy+0.3*dy],[0,0],"-",alpha=1,color="black",lw=1)
        plt.text(self.b/2+self.z2*y+dy+0.3*dy,0.5*y,"y="+str(round(y,3))+"  m",style="italic")
            
        plt.plot([-self.b/2-self.z1*y,self.b/2+self.z2*y],[y+dy,y+dy],"-",alpha=1,color="black",lw=1)
        plt.plot([-self.b/2-self.z1*y,-self.b/2-self.z1*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.plot([self.b/2+self.z2*y,self.b/2+self.z2*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.text(0,y+2*dy,str(round(T,3))+"  m",style="italic")
        plt.text(0-0.1*T,y+2*dy,"T=",style="italic")

        if (self.b==0):
            plt.title("Tirante normal sección triangular",style="italic",fontweight="bold")
        elif (self.z1==0 and self.z2==0):
            plt.title("Tirante normal sección rectangular",style="italic",fontweight="bold")
        elif (self.z1!=0 or self.z2!=0):
            plt.title("Tirante normal sección trapezoidal",style="italic",fontweight="bold") 
        plt.show()

        
#res=Tiranten_ty(z1=1,z2=1,n=0.014,sf=0.001,Q=0.5,b=0.6)
#res.results()
#res.show_section()
#res.show_section3d()


class Tiranten_cy:
    def __init__(self,Q,D,n,sf):
        self.Q=Q
        self.D=D
        self.n=n
        self.sf=sf
    def results(self,decimal_number=6,tol=0.0001,show_iteration="FALSE",show_results="TRUE"):
        error=50
        N=50
        lis_QY=lisQY(self.D,self.n,self.sf,N)
        QQ=lis_QY[0]
        Y=lis_QY[1]
        Qmax=max(QQ)
        Ymax=Y[QQ.index(Qmax)]
        lis_Y=[self.D/2,Ymax]
        lis_Cau=[]
        for i in range(2):
            lis_Cau.append(caudal_sc(lis_Y[i],self.D,self.sf,self.n))
        if (0<self.Q and self.Q<lis_Cau[0]):
            y0=self.D/2
        if (lis_Cau[0]<self.Q and  self.Q<lis_Cau[1]):
            y0=self.D/2+0.5*(Ymax-self.D/2)
        if (show_iteration=="TRUE"):    
            print("yi[m]".center(10),"A[m2]".center(10),"dA".center(10),"P".center(10),"dP".center(10),
                  "f(y)".center(10),"df(y)".center(10),"yi+1[m]".center(10),"error".center(10))
            print(" ==================================================================================================")
        lis_y0=[]
        lis_A=[]
        lis_dA=[]
        lis_P=[]
        lis_dP=[]
        lis_fy=[]
        lis_dfy=[]
        lis_y1=[]
        lis_error=[]
        while (error>tol):
            A=0.25*self.D**2*acos(1-2*y0/self.D)-0.5*self.D**2*(y0/self.D-y0**2/self.D**2)**0.5*(1-2*y0/self.D)
            P=self.D*acos(1-2*y0/self.D)
            dA=0.25*self.D**2*((2/self.D)/(1-(1-2*y0/self.D)**2)**0.5)-0.5*self.D**2*((1-2*y0/self.D)*((1/self.D-2*y0/self.D**2)/(2*(y0/self.D-y0**2/self.D**2)**0.5))+(y0/self.D-y0**2/self.D**2)**0.5*(-2/self.D))
            dP=self.D*(2/self.D)/(1-(1-2*y0/self.D)**2)**0.5
            fy=A**(5/3)/P**(2/3)-(self.Q*self.n)/(self.sf**(0.5))
            dfy=(5/3)*A**(2/3)*P**(-2/3)*dA-(2/3)*A**(5/3)*P**(-5/3)*dP
            y1=y0-fy/dfy
            error=abs(y1-y0)
            lis_y0.append(y0)
            lis_A.append(A)
            lis_dA.append(dA)
            lis_P.append(P)
            lis_dP.append(dP)
            lis_fy.append(fy)
            lis_dfy.append(dfy)
            lis_y1.append(y1)
            lis_error.append(error)
            if (show_iteration=="TRUE"): 
                print(str(round(y0,decimal_number)).center(10),str(round(A,decimal_number)).center(10),
                str(round(dA,decimal_number)).center(10),str(round(P,decimal_number)).center(10),str(round(dP,decimal_number)).center(10),
                str(round(fy,decimal_number)).center(10),str(round(dfy,decimal_number)).center(10),str(round(y1,decimal_number)).center(10),
                str(round(error,decimal_number)).center(10))

            y0=y1
        if (show_iteration=="TRUE"): 
            print(" ==================================================================================================")    
        g=9.807
        ff=len(lis_y0)
        mat_ite=makemat(0,ff,9)
        for i in range(ff):
            mat_ite[i][0]=lis_y0
            mat_ite[i][1]=lis_A
            mat_ite[i][2]=lis_dA
            mat_ite[i][3]=lis_P
            mat_ite[i][4]=lis_dP
            mat_ite[i][5]=lis_fy
            mat_ite[i][6]=lis_dfy
            mat_ite[i][7]=lis_y1
            mat_ite[i][8]=lis_error
        
        teta=2*acos(1-(2*y1)/self.D)
        T=self.D*sin(teta/2)
        F=(self.Q/A)/(9.807*A/(abs(T)))**0.5
        E= y1+(self.Q**2)/(2*g*A**2)
        V= self.Q/A
        R=A/P
        if (show_results=="TRUE"):
            print(" Tirante normal          yn[m] :",y1)
            print(" Área hidráulica         A[m2] :",A)
            print(" Espejo de agua           T[m] :",T)
            print(" Número de Froude            F :",F)
            print(" Perímetro mojado         P[m] :",P)
            print(" Radio hidráulico         R[m] :",R)
            print(" Velocidad              V[m/s] :",V)
            print(" Energía específica E[m-kg/kg] :",E)
        return [y1,A,T,F,P,R,V,E,mat_ite]
        
    def show_section(self,outline_color="silver",fill_color="cyan"):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        thetar=2*acos(1-(2*y)/(self.D))
        thetag=(180/pi)*thetar
        theta11=180+(90-thetag/2)
        theta22=270+(thetag/2)
        fig1 = plt.figure()
    ####################
        ax1 = fig1.add_subplot(111, aspect='equal',xticks=[], yticks=[])
    ######################
        ax1.add_patch(patches.Wedge(center=(0,0), r=self.D/2+self.D/35, theta1=0
                                    , theta2=360, facecolor=outline_color, label="Test"))
        ax1.add_patch(patches.Wedge(center=(0,0), r=self.D/2, theta1=0
                                   , theta2=360, facecolor="white", label="Test"))
#######################
        ax1.add_patch(patches.Wedge(center=(0,0), r=self.D/2, theta1=theta11
                                   , theta2=theta22, facecolor=fill_color, label="Test"))
        plt.xlim(-self.D/2-self.D/2, self.D/2+self.D/2)
        plt.ylim(-self.D/2-self.D/3, self.D/2+self.D/3)
        array=np.array([[0,0],[-self.D/2*sin(thetar/2),-(self.D/2-y)],[self.D/2*sin(thetar/2),-(self.D/2-y)]])
        if (y<=self.D/2):
            shape=Polygon(array,color="white")
        elif (y>self.D/2):
            shape=Polygon(array,color=fill_color)
        ax1.add_patch(shape)
        plt.title("Tirante normal sección circular",style="italic",fontweight="bold")
    
        plt.plot([self.D/2+0.07*self.D,self.D/2+0.07*self.D],[-self.D/2,-self.D/2+y], "-",alpha=1,color="black",lw=1)
        plt.plot([self.D/2+0.07*self.D-0.02*self.D,self.D/2+0.07*self.D+0.02*self.D],[-self.D/2+y,-self.D/2+y], "-",alpha=1,color="black",lw=1)
        plt.plot([self.D/2+0.07*self.D-0.02*self.D,self.D/2+0.07*self.D+0.02*self.D],[-self.D/2,-self.D/2], "-",alpha=1,color="black",lw=1)
        plt.text(self.D/2+0.07*self.D+0.03*self.D,-self.D/2+0.5*y,"y="+str(round(y,3))+" m",style="italic")

        plt.plot([-self.D/2-0.07*self.D,-self.D/2-0.07*self.D],[-self.D/2,self.D/2], "-",alpha=1,color="black",lw=1)
        plt.plot([-self.D/2-0.07*self.D-0.02*self.D,-self.D/2-0.07*self.D+0.02*self.D],[self.D/2,self.D/2], "-",alpha=1,color="black",lw=1)
        plt.plot([-self.D/2-0.07*self.D-0.02*self.D,-self.D/2-0.07*self.D+0.02*self.D],[-self.D/2,-self.D/2], "-",alpha=1,color="black",lw=1)
        plt.text(-self.D/2-0.07*self.D-0.4*self.D,-self.D/2+0.5*self.D,"D="+str(round(self.D,3))+" m",style="italic")

        teta=2*acos(1-(2*y)/self.D)
        T=self.D*sin(teta/2)
        plt.plot([-0.5*T,0.5*T],[self.D/2+0.07*self.D,self.D/2+0.07*self.D], "-",alpha=1,color="black",lw=1)
        plt.plot([-0.5*T,-0.5*T],[self.D/2+0.07*self.D-0.02*self.D,self.D/2+0.07*self.D+0.02*self.D], "-",alpha=1,color="black",lw=1)
        plt.plot([0.5*T,0.5*T],[self.D/2+0.07*self.D-0.02*self.D,self.D/2+0.07*self.D+0.02*self.D], "-",alpha=1,color="black",lw=1)
        plt.text(0,self.D/2+0.07*self.D+0.03*self.D,str(round(T,3))+" m",style="italic")
        plt.text(-0.18*self.D,self.D/2+0.07*self.D+0.03*self.D,"T=",style="italic")    
        N=50
        lis_QY=lisQY(self.D,self.n,self.sf,N)
        lis_VY=lisVY(self.D,self.n,self.sf,N)
        QQ=lis_QY[0]
        VV=lis_VY[0]     
        Y=lis_QY[1]
        Qmax=max(QQ)
        Vmax=max(VV)
        plt.text(-self.D/2-self.D/2.5,-self.D/2-self.D/5,"Qmax="+str(round(Qmax,3))+" m3/s",style="italic")
        plt.text(-self.D/2-self.D/2.5,-self.D/2-self.D/3.5,"Vmax="+str(round(Vmax,3))+" m/s",style="italic")
        plt.show()
#rr=Tiranten_cy(Q=0.1,D=0.6,n=0.014,sf=0.001)
#rr.results(show_iteration="TRUE")
#rr.show_section(outline_color="red",fill_color="orange")
class Max_efficiency:
    def __init__(self,Q,z,n,sf):
        self.Q=Q
        self.z=z
        self.n=n
        self.sf=sf
    def results(self,y0=0.5,decimal_number=6,tol=0.0001,show_iteration="FALSE",show_results="TRUE"):
        error=500
        if (show_iteration=="TRUE"):
            print("yi[m]".center(10),"b[m]".center(10),"A[m2]".center(10),
                  "dA".center(10),"P".center(10),"dP".center(10),
                  "f(y)".center(10),"df(y)".center(10),"yi+1[m]".center(10),
                  "error".center(10))
            print(" ===========================================================================================================")
        lis_y0=[]
        lis_b=[]
        lis_A=[]
        lis_dA=[]
        lis_P=[]
        lis_dP=[]
        lis_fy=[]
        lis_dfy=[]
        lis_y1=[]
        lis_error=[]
        while (error>tol):
            b=2*((1+self.z**2)**0.5-self.z)*y0
            A=b*y0+self.z*y0**2
            P=b+2*y0*(self.z**2+1)**0.5
            db=2*((1+self.z)**0.5-self.z)
            dA=y0*db+b+2*self.z*y0
            dP=db+2*(self.z**2+1)**0.5
            fy=(A**(5/3))/(P**(2/3))-(self.Q*self.n)/(self.sf**0.5)
            dfy=A**(5/3)*(-2/3*P**(-5/3)*dP)+P**(-2/3)*(5/3*A**(2/3)*dA)
            y1=y0-fy/dfy
            error=abs(y1-y0)
            lis_y0.append(y0)
            lis_b.append(b)
            lis_A.append(A)
            lis_dA.append(dA)
            lis_P.append(P)
            lis_dP.append(dP)
            lis_fy.append(fy)
            lis_dfy.append(dfy)
            lis_y1.append(y1)
            lis_error.append(error)
            if (show_iteration=="TRUE"):
                print(str(round(y0,decimal_number)).center(10),str(round(b,decimal_number)).center(10),
                      str(round(A,decimal_number)).center(10),
                      str(round(dA,decimal_number)).center(10),str(round(P,decimal_number)).center(10),str(round(dP,decimal_number)).center(10),
                      str(round(fy,decimal_number)).center(10),str(round(dfy,decimal_number)).center(10),str(round(y1,decimal_number)).center(10),
                      str(round(error,decimal_number)).center(10))
            y0=y1
        if (show_iteration=="TRUE"):
            print(" ===========================================================================================================")

        ff=len(lis_y0)
        mat_ite=makemat(0,ff,10)
        for i in range(ff):
            mat_ite[i][0]=lis_y0
            mat_ite[i][1]=lis_b
            mat_ite[i][2]=lis_A
            mat_ite[i][3]=lis_dA
            mat_ite[i][4]=lis_P
            mat_ite[i][5]=lis_dP
            mat_ite[i][6]=lis_fy
            mat_ite[i][7]=lis_dfy
            mat_ite[i][8]=lis_y1
            mat_ite[i][9]=lis_error


        g=9.807
        T=self.z*y0+b+self.z*y0
        R=A/P
        V=self.Q/A
        F=(V)/((g*(A)/(T))**0.5)
        E=y0+(V**2)/(2*g)
        if (show_results=="TRUE"):
            print(" Tirante normal          yn[m] :",y0)
            print(" Ancho de solera          b[m] :",b)
            print(" Área hidráulica         A[m2] :",A)
            print(" Espejo de agua           T[m] :",T)
            print(" Número de Froude            F :",F)
            print(" Perímetro mojado         P[m] :",P)
            print(" Radio hidráulico         R[m] :",R)
            print(" Velocidad              V[m/s] :",V)
            print(" Energía específica E[m-kg/kg] :",E)
        return [y1,b,A,T,F,P,R,V,E,mat_ite]
    def show_section(self,outline_color="silver",fill_color="cyan",width_size=7,height_size=5):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        b=res[1]
        z1=self.z
        z2=self.z
        T=z1*y+b+z2*y
        fig1=plt.figure(figsize=(width_size,height_size)) 
        ax1=fig1.add_subplot(xticks=[], yticks=[])
        dx1=-b/2-1.25*y*z1 
        dx2=b/2+1.25*y*z2
        plt.xlim(dx1-0.25*T,dx2+0.25*T)
        plt.ylim(-0.7*T+0.5*y,0.7*T+0.5*y)
        yyy=0
        xxx=0
        array=np.array([[-b/2+xxx,0+yyy],[-b/2-z1*y+xxx,y+yyy],
                          [b/2+z2*y+xxx,y+yyy],[b/2+xxx,0+yyy]])
        shape=Polygon(array,color=fill_color)
        ax1.add_patch(shape)
        plt.plot([-b/2-1.25*y*z1+xxx,-b/2+xxx,b/2+xxx,b/2+1.25*y*z2+xxx],
                   [1.25*y+yyy,0+yyy,0+yyy,1.25*y+yyy],"-",color=outline_color,lw=6)
        dy=0.05*abs(-dx1+dx2)
        if (b!=0):
            plt.plot([-b/2,b/2],[-dy,-dy],"-",alpha=1,color="black",lw=1)
            plt.plot([-b/2,-b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.plot([b/2,b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.text(0,-3*dy,str(round(b,3))+"  m",style="italic")
            plt.text(0-0.25*b,-3*dy,"b=",style="italic")

        plt.plot([b/2+z2*y+dy,b/2+z2*y+dy],[0,y],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y+dy-0.3*dy ,b/2+z2*y+dy+0.3*dy],[y,y],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y+dy-0.3*dy ,b/2+z2*y+dy+0.3*dy],[0,0],"-",alpha=1,color="black",lw=1)
        plt.text(b/2+z2*y+dy+0.3*dy,0.5*y,"y="+str(round(y,3))+"  m",style="italic")
            
        plt.plot([-b/2-z1*y,b/2+z2*y],[y+dy,y+dy],"-",alpha=1,color="black",lw=1)
        plt.plot([-b/2-z1*y,-b/2-z1*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y,b/2+z2*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.text(0,y+2*dy,str(round(T,3))+"  m",style="italic")
        plt.text(0-0.1*T,y+2*dy,"T=",style="italic")

        #if (b==0):
        plt.title("Sección de máxima eficiencia hidráulica",style="italic",fontweight="bold")
        #elif (z1==0 and z2==0):
        #    plt.title("Tirante normal sección rectangular",style="italic",fontweight="bold")
        #elif (z1!=0 or z2!=0):
        #    plt.title("Tirante normal sección trapezoidal",style="italic",fontweight="bold") 
        plt.show()
    def show_section3d(self,lis_col=["red","lightcoral","cyan","blue"]):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        b=res[1]
        z1=self.z
        z2=self.z
        L=10
        rr=0.25
        h=y+rr*y
        points = np.array([[-0.5*L,-0.5*b,0],
                          [0.5*L, -0.5*b, 0 ],
                          [0.5*L,-0.5*b-z1*h, h],
                          [-0.5*L,-0.5*b-z1*h,h],
                          [-0.5*L,0.5*b,0],
                          [0.5*L, 0.5*b, 0 ],
                          [0.5*L,0.5*b+z2*h, h],
                          [-0.5*L,0.5*b+z2*h,h] 
                           ])
        #ancho=(z1*h+b+z2*h)/15
        ancho=0.2*(z1*y+b+z2*y)
        anchov=0.2*(y)
        points2 = np.array([[-0.5*L,-0.5*b-ancho-z1*ancho,-anchov],
                          [0.5*L, -0.5*b-ancho-z1*ancho, -anchov],
                          [0.5*L,-0.5*b-z1*h-ancho, h],
                          [-0.5*L,-0.5*b-z1*h-ancho,h],
                          [-0.5*L,0.5*b+ancho+z2*ancho,-anchov],
                          [0.5*L, 0.5*b+ancho+z2*ancho, -anchov],
                          [0.5*L,0.5*b+z2*h+ancho, h],
                          [-0.5*L,0.5*b+z2*h+ancho,h] 
                           ])
    
        points1 = np.array([[-0.5*L,-0.5*b,0],
                          [0.5*L, -0.5*b, 0 ],
                          [0.5*L,-0.5*b-z1*y, y],
                          [-0.5*L,-0.5*b-z1*y,y],
                          [-0.5*L,0.5*b,0],
                          [0.5*L, 0.5*b, 0 ],
                          [0.5*L,0.5*b+z2*y, y],
                          [-0.5*L,0.5*b+z2*y,y] 
                           ])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        r = [-5,5]
        X, Y = np.meshgrid(r,r)
        ax.set(xlim3d=(-0.5*L,0.5*L), xlabel='L')
        ax.set(ylim3d=(-0.5*b-z1*h-ancho-0.1*h,0.5*b+z2*h+ancho+0.1*h ),
               ylabel='')
        ax.set(zlim3d=(-ancho-0.3*h,h+0.3*h), zlabel='y')
        # vértices 
        #ax.scatter3D(points[:,0],points[:,1],points[:,2]) 
        #Lista de polígonos de los lados de la figura
        verts=[[points[0],points[1],points[2],points[3]],
               [points[4],points[5],points[6],points[7]],
               [points[0],points[1],points[5],points[4]]
               ]
        verts2=[[points2[0],points2[1],points2[2],points2[3]],
               [points2[4],points2[5],points2[6],points2[7]],
               [points2[0],points2[1],points2[5],points2[4]],
               [points[2],points[3],points2[3],points2[2]],
               [points[6],points[7],points2[7],points2[6]]  
               ]
        verts1=[[points1[0],points1[1],points1[2],points1[3]],
                [points1[4],points1[5],points1[6],points1[7]],
                [points1[0],points1[1],points1[5],points1[4]],
                [points1[0],points1[4],points1[7],points1[3]],
                [points1[1],points1[5],points1[6],points1[2]],
                [points1[2],points1[6],points1[7],points1[3]]
                ]
    #lados de la trama
        ax.add_collection3d(Poly3DCollection(verts,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))

        ax.add_collection3d(Poly3DCollection(verts2,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))
        ax.add_collection3d(Poly3DCollection(verts1,
                                             facecolors=lis_col[2],linewidths=1,
                                             edgecolors=lis_col[3], alpha=0.5))
        #if (self.b==0):
        #    plt.title("Tirante normal sección triangular")
        #elif (self.z1==0 and self.z2==0):
        plt.title("Sección de máxima eficiencia hidráulica")
        #elif (self.z1!=0 or self.z2!=0):
        #    plt.title("Tirante normal sección trapezoidal")    
        plt.show()        

#rr=Max_efficiency(sf=0.001,Q=0.6,z=0.5,n=0.014)
#rr.results(show_iteration="TRUE")
#rr.show_section3d()
class Min_infiltration:
    def __init__(self,Q,z,n,sf):
        self.Q=Q
        self.z=z
        self.n=n
        self.sf=sf
    def results(self,y0=0.5,decimal_number=6,tol=0.0001,show_iteration="FALSE",show_results="TRUE"):
        error=500
        if (show_iteration=="TRUE"):
            print("yi[m]".center(10),"b[m]".center(10),"A[m2]".center(10),
                  "dA".center(10),"P".center(10),"dP".center(10),
                  "f(y)".center(10),"df(y)".center(10),"yi+1[m]".center(10),
                  "error".center(10))
            print(" ===========================================================================================================")
        lis_y0=[]
        lis_b=[]
        lis_A=[]
        lis_dA=[]
        lis_P=[]
        lis_dP=[]
        lis_fy=[]
        lis_dfy=[]
        lis_y1=[]
        lis_error=[]
        while (error>tol):
            b=4*((1+self.z**2)**0.5-self.z)*y0
            A=b*y0+self.z*y0**2
            P=b+2*y0*(self.z**2+1)**0.5
            db=2*((1+self.z)**0.5-self.z)
            dA=y0*db+b+2*self.z*y0
            dP=db+2*(self.z**2+1)**0.5
            fy=(A**(5/3))/(P**(2/3))-(self.Q*self.n)/(self.sf**0.5)
            dfy=A**(5/3)*(-2/3*P**(-5/3)*dP)+P**(-2/3)*(5/3*A**(2/3)*dA)
            y1=y0-fy/dfy
            error=abs(y1-y0)
            lis_y0.append(y0)
            lis_b.append(b)
            lis_A.append(A)
            lis_dA.append(dA)
            lis_P.append(P)
            lis_dP.append(dP)
            lis_fy.append(fy)
            lis_dfy.append(dfy)
            lis_y1.append(y1)
            lis_error.append(error)
            if (show_iteration=="TRUE"):
                print(str(round(y0,decimal_number)).center(10),str(round(b,decimal_number)).center(10),
                      str(round(A,decimal_number)).center(10),
                      str(round(dA,decimal_number)).center(10),str(round(P,decimal_number)).center(10),str(round(dP,decimal_number)).center(10),
                      str(round(fy,decimal_number)).center(10),str(round(dfy,decimal_number)).center(10),str(round(y1,decimal_number)).center(10),
                      str(round(error,decimal_number)).center(10))
            y0=y1
        if (show_iteration=="TRUE"):
            print(" ===========================================================================================================")

        ff=len(lis_y0)
        mat_ite=makemat(0,ff,10)
        for i in range(ff):
            mat_ite[i][0]=lis_y0
            mat_ite[i][1]=lis_b
            mat_ite[i][2]=lis_A
            mat_ite[i][3]=lis_dA
            mat_ite[i][4]=lis_P
            mat_ite[i][5]=lis_dP
            mat_ite[i][6]=lis_fy
            mat_ite[i][7]=lis_dfy
            mat_ite[i][8]=lis_y1
            mat_ite[i][9]=lis_error


        g=9.807
        T=self.z*y0+b+self.z*y0
        R=A/P
        V=self.Q/A
        F=(V)/((g*(A)/(T))**0.5)
        E=y0+(V**2)/(2*g)
        if (show_results=="TRUE"):
            print(" Tirante normal          yn[m] :",y0)
            print(" Ancho de solera          b[m] :",b)
            print(" Área hidráulica         A[m2] :",A)
            print(" Espejo de agua           T[m] :",T)
            print(" Número de Froude            F :",F)
            print(" Perímetro mojado         P[m] :",P)
            print(" Radio hidráulico         R[m] :",R)
            print(" Velocidad              V[m/s] :",V)
            print(" Energía específica E[m-kg/kg] :",E)
        return [y1,b,A,T,F,P,R,V,E,mat_ite]
    def show_section(self,outline_color="silver",fill_color="cyan",width_size=7,height_size=5):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        b=res[1]
        z1=self.z
        z2=self.z
        T=z1*y+b+z2*y
        fig1=plt.figure(figsize=(width_size,height_size)) 
        ax1=fig1.add_subplot(xticks=[], yticks=[])
        dx1=-b/2-1.25*y*z1 
        dx2=b/2+1.25*y*z2
        plt.xlim(dx1-0.25*T,dx2+0.25*T)
        plt.ylim(-0.7*T+0.5*y,0.7*T+0.5*y)
        yyy=0
        xxx=0
        array=np.array([[-b/2+xxx,0+yyy],[-b/2-z1*y+xxx,y+yyy],
                          [b/2+z2*y+xxx,y+yyy],[b/2+xxx,0+yyy]])
        shape=Polygon(array,color=fill_color)
        ax1.add_patch(shape)
        plt.plot([-b/2-1.25*y*z1+xxx,-b/2+xxx,b/2+xxx,b/2+1.25*y*z2+xxx],
                   [1.25*y+yyy,0+yyy,0+yyy,1.25*y+yyy],"-",color=outline_color,lw=6)
        dy=0.05*abs(-dx1+dx2)
        if (b!=0):
            plt.plot([-b/2,b/2],[-dy,-dy],"-",alpha=1,color="black",lw=1)
            plt.plot([-b/2,-b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.plot([b/2,b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.text(0,-3*dy,str(round(b,3))+"  m",style="italic")
            plt.text(0-0.25*b,-3*dy,"b=",style="italic")

        plt.plot([b/2+z2*y+dy,b/2+z2*y+dy],[0,y],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y+dy-0.3*dy ,b/2+z2*y+dy+0.3*dy],[y,y],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y+dy-0.3*dy ,b/2+z2*y+dy+0.3*dy],[0,0],"-",alpha=1,color="black",lw=1)
        plt.text(b/2+z2*y+dy+0.3*dy,0.5*y,"y="+str(round(y,3))+"  m",style="italic")
            
        plt.plot([-b/2-z1*y,b/2+z2*y],[y+dy,y+dy],"-",alpha=1,color="black",lw=1)
        plt.plot([-b/2-z1*y,-b/2-z1*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y,b/2+z2*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.text(0,y+2*dy,str(round(T,3))+"  m",style="italic")
        plt.text(0-0.1*T,y+2*dy,"T=",style="italic")

        #if (b==0):
        plt.title("Sección de mínima infiltración",style="italic",fontweight="bold")
        #elif (z1==0 and z2==0):
        #    plt.title("Tirante normal sección rectangular",style="italic",fontweight="bold")
        #elif (z1!=0 or z2!=0):
        #    plt.title("Tirante normal sección trapezoidal",style="italic",fontweight="bold") 
        plt.show()
    def show_section3d(self,lis_col=["red","lightcoral","cyan","blue"]):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        b=res[1]
        z1=self.z
        z2=self.z
        L=10
        rr=0.25
        h=y+rr*y
        points = np.array([[-0.5*L,-0.5*b,0],
                          [0.5*L, -0.5*b, 0 ],
                          [0.5*L,-0.5*b-z1*h, h],
                          [-0.5*L,-0.5*b-z1*h,h],
                          [-0.5*L,0.5*b,0],
                          [0.5*L, 0.5*b, 0 ],
                          [0.5*L,0.5*b+z2*h, h],
                          [-0.5*L,0.5*b+z2*h,h] 
                           ])
        #ancho=(z1*h+b+z2*h)/15
        ancho=0.2*(z1*y+b+z2*y)
        anchov=0.2*(y)
        points2 = np.array([[-0.5*L,-0.5*b-ancho-z1*ancho,-anchov],
                          [0.5*L, -0.5*b-ancho-z1*ancho, -anchov],
                          [0.5*L,-0.5*b-z1*h-ancho, h],
                          [-0.5*L,-0.5*b-z1*h-ancho,h],
                          [-0.5*L,0.5*b+ancho+z2*ancho,-anchov],
                          [0.5*L, 0.5*b+ancho+z2*ancho, -anchov],
                          [0.5*L,0.5*b+z2*h+ancho, h],
                          [-0.5*L,0.5*b+z2*h+ancho,h] 
                           ])
    
        points1 = np.array([[-0.5*L,-0.5*b,0],
                          [0.5*L, -0.5*b, 0 ],
                          [0.5*L,-0.5*b-z1*y, y],
                          [-0.5*L,-0.5*b-z1*y,y],
                          [-0.5*L,0.5*b,0],
                          [0.5*L, 0.5*b, 0 ],
                          [0.5*L,0.5*b+z2*y, y],
                          [-0.5*L,0.5*b+z2*y,y] 
                           ])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        r = [-5,5]
        X, Y = np.meshgrid(r,r)
        ax.set(xlim3d=(-0.5*L,0.5*L), xlabel='L')
        ax.set(ylim3d=(-0.5*b-z1*h-ancho-0.1*h,0.5*b+z2*h+ancho+0.1*h ),
               ylabel='')
        ax.set(zlim3d=(-ancho-0.3*h,h+0.3*h), zlabel='y')
        # vértices 
        #ax.scatter3D(points[:,0],points[:,1],points[:,2]) 
        #Lista de polígonos de los lados de la figura
        verts=[[points[0],points[1],points[2],points[3]],
               [points[4],points[5],points[6],points[7]],
               [points[0],points[1],points[5],points[4]]
               ]
        verts2=[[points2[0],points2[1],points2[2],points2[3]],
               [points2[4],points2[5],points2[6],points2[7]],
               [points2[0],points2[1],points2[5],points2[4]],
               [points[2],points[3],points2[3],points2[2]],
               [points[6],points[7],points2[7],points2[6]]  
               ]
        verts1=[[points1[0],points1[1],points1[2],points1[3]],
                [points1[4],points1[5],points1[6],points1[7]],
                [points1[0],points1[1],points1[5],points1[4]],
                [points1[0],points1[4],points1[7],points1[3]],
                [points1[1],points1[5],points1[6],points1[2]],
                [points1[2],points1[6],points1[7],points1[3]]
                ]
    #lados de la trama
        ax.add_collection3d(Poly3DCollection(verts,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))

        ax.add_collection3d(Poly3DCollection(verts2,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))
        ax.add_collection3d(Poly3DCollection(verts1,
                                             facecolors=lis_col[2],linewidths=1,
                                             edgecolors=lis_col[3], alpha=0.5))
        #if (self.b==0):
        #    plt.title("Tirante normal sección triangular")
        #elif (self.z1==0 and self.z2==0):
        plt.title("Sección de mínima infiltración")
        #elif (self.z1!=0 or self.z2!=0):
        #    plt.title("Tirante normal sección trapezoidal")    
        plt.show()
#rr=Min_infiltration(Q=0.6,z=0.5,n=0.014,sf=0.001)
#rr.results(show_iteration="TRUE")
#rr.show_section()
class RP_HortonEinstein:
    def __init__(self,mat_dat,sf,h0,n0,nf):
        self.mat_dat=mat_dat
        self.sf=sf
        self.h0=h0
        self.n0=n0
        self.nf=nf
    def results(self,show_results="TRUE"):    
        lis_x=[0,0]
        lis_y=[0,-self.h0]
        xx=0
        for i in range(len(self.mat_dat)):
            xx=xx+self.mat_dat[i][0]
            lis_x.append(xx)
            lis_y.append(-self.mat_dat[i][1])
        lis_x.append(xx)
        lis_y.append(0)
        lis_h=[self.h0]
        for i in range(len(self.mat_dat)):
            lis_h.append(self.mat_dat[i][1])
        nn=len(lis_h)
        lis_a=[]
        lis_p=[]
        for i in range(nn-1):
            lis_a.append(0.5*(lis_h[i]+lis_h[i+1])*self.mat_dat[i][0])
            lis_p.append(((lis_h[i+1]-lis_h[i])**2+self.mat_dat[i][0]**2)**0.5)
        P=sum(lis_p)+self.mat_dat[nn-2][1]+self.h0
        sumad=0
        for i in range(len(self.mat_dat)):
            sumad=sumad+lis_p[i]*self.mat_dat[i][2]**1.5   
        sumad=sumad+(self.h0*self.n0**1.5+self.mat_dat[nn-2][1]*self.nf**1.5)
        np=(sumad**(2/3))/(P**(2/3))
        A=sum(lis_a)
        Q=(A**(5/3)*self.sf**(1/2))/(P**(2/3)*np)
        V=Q/A
        T=max(lis_x)
        g=9.807
        F=(V)/((g*((A)/(T)))**0.5)
        if (F<1):
            tex="Subcrítico"
        elif (F==1):
            tex="Crítico"
        elif (F>1):
            tex="Supercrítico"
        if (show_results=="TRUE"):    
            print("==========Rugosidad ponderada Horton y Einstein =========")
            print("Caudal                 (Q)[m3/s] :",Q)
            print("Área hidráulica          (A)[m2] :",A)
            print("Perímetro mojado          (P)[m] :",P)
            print("Rugosidad ponderada         (np) :",np)
            print("Velocidad               (V)[m/s] :",V)
            print("Espejo de agua            (T)[m] :",T)
            print("Número de Froude             (F) :",F)
            print("Tipo de flujo                    :",tex)
            print("=========================================================")
        return [Q,A,P,np,V,T,F,lis_x,lis_y]
    #def plot(lis_x,lis_y,lis_col,tex):
    def show_section(self,lis_col=["cyan","red","magenta"],alpha=0.9):
        res=self.results(show_results="FALSE")
        lis_x=res[7]
        lis_y=res[8]
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        xx=max(lis_x)
        plt.xlim(-0.05*xx,xx+0.05*xx)
        limy=min(lis_y)
        plt.ylim(+limy+0.1*limy,0-0.1*limy)
        datxy=[]
        for i in range(len(lis_x)):
            datxy.append([lis_x[i],lis_y[i]])
        shape=Polygon(datxy,color=lis_col[0],alpha=alpha)
        plt.plot(lis_x,lis_y,"-",color=lis_col[1])
        plt.plot(lis_x,lis_y,"o",color=lis_col[2])
        ax1.add_patch(shape)
        plt.grid()
        plt.title("Rugosidad ponderada Horton y Einstein")
        plt.show()



#mat_dat=[[1.5,1.4,0.014],
#         [3,0.5,0.015]]
#sf=0.0005 #Pendiente de la linea de energia
#h0=1.4 #Profundidad en el margen (h0)[m]
#n0=0.016 #Coeficiente de Manning en el margen (n0)
#nf=0.016
#rr=RP_HortonEinstein(mat_dat=mat_dat,sf=sf,h0=h0,n0=n0,nf=nf)
#rr.results()
#rr.show_section(lis_col=["blue","red","yellow"],alpha=0.1)
class RP_Lotter:
    def __init__(self,mat_dat,sf,h0,n0,nf):
        self.mat_dat=mat_dat
        self.sf=sf
        self.h0=h0
        self.n0=n0
        self.nf=nf
    def results(self,show_results="TRUE"):    
        lis_x=[0,0]
        lis_y=[0,-self.h0]
        xx=0
        for i in range(len(self.mat_dat)):
            xx=xx+self.mat_dat[i][0]
            lis_x.append(xx)
            lis_y.append(-self.mat_dat[i][1])
        lis_x.append(xx)
        lis_y.append(0)
        lis_h=[self.h0]
        for i in range(len(self.mat_dat)):
            lis_h.append(self.mat_dat[i][1])
        nn=len(lis_h)
        lis_a=[]
        lis_p=[]
        lis_R=[]
        for i in range(nn-1):
            lis_a.append(0.5*(lis_h[i]+lis_h[i+1])*self.mat_dat[i][0])
            lis_p.append(((lis_h[i+1]-lis_h[i])**2+self.mat_dat[i][0]**2)**0.5)
            lis_R.append(lis_a[i]/lis_p[i])
        P=sum(lis_p)+self.mat_dat[nn-2][1]+h0
        R=sum(lis_R)
        sumad=0
        for i in range(len(self.mat_dat)):
            sumad=sumad+(lis_p[i]*lis_R[i]**(5/3))/(self.mat_dat[i][2])   
        np=(P*R**(5/3))/(sumad)
        A=sum(lis_a)
        Q=(A**(5/3)*self.sf**(1/2))/(P**(2/3)*np)
        V=Q/A
        T=max(lis_x)
        g=9.807
        F=(V)/((g*((A)/(T)))**0.5)
        if (F<1):
            tex="Subcrítico"
        elif (F==1):
            tex="Crítico"
        elif (F>1):
            tex="Supercrítico"
        print("================Rugosidad ponderada Lotter ==============")
        print("Caudal                 (Q)[m3/s] :",Q)
        print("Área hidráulica          (A)[m2] :",A)
        print("Perímetro mojado          (P)[m] :",P)
        print("Rugosidad ponderada         (np) :",np)
        print("Velocidad               (V)[m/s] :",V)
        print("Espejo de agua            (T)[m] :",T)
        print("Número de Froude             (F) :",F)
        print("Tipo de flujo                    :",tex)
        print("=========================================================")
        return [Q,A,P,np,V,T,F,lis_x,lis_y]
    #def plot(lis_x,lis_y,lis_col,tex):
    def show_section(self,lis_col=["cyan","red","magenta"],alpha=0.9):
        res=self.results(show_results="FALSE")
        lis_x=res[7]
        lis_y=res[8]
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        xx=max(lis_x)
        plt.xlim(-0.05*xx,xx+0.05*xx)
        limy=min(lis_y)
        plt.ylim(+limy+0.1*limy,0-0.1*limy)
        datxy=[]
        for i in range(len(lis_x)):
            datxy.append([lis_x[i],lis_y[i]])
        shape=Polygon(datxy,color=lis_col[0],alpha=alpha)
        plt.plot(lis_x,lis_y,"-",color=lis_col[1])
        plt.plot(lis_x,lis_y,"o",color=lis_col[2])
        ax1.add_patch(shape)
        plt.grid()
        plt.title("Rugosidad ponderada Lotter")
        plt.show()

#mat_dat=[[1.5,1.4,0.014],
#         [3,0.5,0.015]]
#sf=0.0005 #Pendiente de la linea de energia
#h0=1.4 #Profundidad en el margen (h0)[m]
#n0=0.016 #Coeficiente de Manning en el margen (n0)
#nf=0.016
#rr=RP_Lotter(mat_dat,sf,h0,n0,nf)
#rr.results()
#rr.show_section(lis_col=["blue","red","yellow"],alpha=0.1)
class RP_Kelkel:
    def __init__(self,mat_dat,sf,h0,n0,nf):
        self.mat_dat=mat_dat
        self.sf=sf
        self.h0=h0
        self.n0=n0
        self.nf=nf
    def results(self,show_results="TRUE"):
        lis_x=[0,0]
        lis_y=[0,-self.h0]
        xx=0
        for i in range(len(self.mat_dat)):
            xx=xx+self.mat_dat[i][0]
            lis_x.append(xx)
            lis_y.append(-self.mat_dat[i][1])
        lis_x.append(xx)
        lis_y.append(0)
        lis_h=[h0]
        for i in range(len(self.mat_dat)):
            lis_h.append(self.mat_dat[i][1])
        nn=len(lis_h)
        lis_a=[]
        lis_p=[]
        for i in range(nn-1):
            lis_a.append(0.5*(lis_h[i]+lis_h[i+1])*self.mat_dat[i][0])
            lis_p.append(((lis_h[i+1]-lis_h[i])**2+self.mat_dat[i][0]**2)**0.5)
        P=sum(lis_p)+self.mat_dat[nn-2][1]+h0
        sumad=0
        for i in range(len(self.mat_dat)):
            sumad=sumad+lis_p[i]/(self.mat_dat[i][2])   
        sumad=sumad+(self.h0/self.n0  +self.mat_dat[nn-2][1]/self.nf)
        np=(P)/(sumad)
        A=sum(lis_a)
        Q=(A**(5/3)*self.sf**(1/2))/(P**(2/3)*np)
        V=Q/A
        T=max(lis_x)
        g=9.807
        F=(V)/((g*((A)/(T)))**0.5)
        if (F<1):
            tex="Subcrítico"
        elif (F==1):
            tex="Crítico"
        elif (F>1):
            tex="Supercrítico"
        print("============== Rugosidad ponderada Kelkel ===============")
        print("Caudal                 (Q)[m3/s] :",Q)
        print("Área hidráulica          (A)[m2] :",A)
        print("Perímetro mojado          (P)[m] :",P)
        print("Rugosidad ponderada         (np) :",np)
        print("Velocidad               (V)[m/s] :",V)
        print("Espejo de agua            (T)[m] :",T)
        print("Número de Froude             (F) :",F)
        print("Tipo de flujo                    :",tex)
        print("=========================================================")
        return [Q,A,P,np,V,T,F,lis_x,lis_y]

    #def plot(lis_x,lis_y,lis_col,tex):
    def show_section(self,lis_col=["cyan","red","magenta"],alpha=0.9):
        res=self.results(show_results="FALSE")
        lis_x=res[7]
        lis_y=res[8]
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        xx=max(lis_x)
        plt.xlim(-0.05*xx,xx+0.05*xx)
        limy=min(lis_y)
        plt.ylim(+limy+0.1*limy,0-0.1*limy)
        datxy=[]
        for i in range(len(lis_x)):
            datxy.append([lis_x[i],lis_y[i]])
        shape=Polygon(datxy,color=lis_col[0],alpha=alpha)
        plt.plot(lis_x,lis_y,"-",color=lis_col[1])
        plt.plot(lis_x,lis_y,"o",color=lis_col[2])
        ax1.add_patch(shape)
        plt.grid()
        plt.title("Rugosidad ponderada Kelkel")
        plt.show()

#mat_dat=[[1.5,1.4,0.014],
#         [3,0.5,0.015]]
#sf=0.0005 #Pendiente de la linea de energia
#h0=1.4 #Profundidad en el margen (h0)[m]
#n0=0.016 #Coeficiente de Manning en el margen (n0)
#nf=0.016
#rr=RP_Kelkel(mat_dat,sf,h0,n0,nf)
#rr.results()
#rr.show_section(lis_col=["blue","red","yellow"],alpha=0.1)
class RP_Yen:
    def __init__(self,mat_dat,sf,h0,n0,nf):
        self.mat_dat=mat_dat
        self.sf=sf
        self.h0=h0
        self.n0=n0
        self.nf=nf
    def results(self,show_results="TRUE"):
        lis_x=[0,0]
        lis_y=[0,-self.h0]
        xx=0
        for i in range(len(self.mat_dat)):
            xx=xx+self.mat_dat[i][0]
            lis_x.append(xx)
            lis_y.append(-self.mat_dat[i][1])
        lis_x.append(xx)
        lis_y.append(0)
        lis_h=[self.h0]
        for i in range(len(self.mat_dat)):
            lis_h.append(self.mat_dat[i][1])
        nn=len(lis_h)
        lis_a=[]
        lis_p=[]
        lis_R=[]
        for i in range(nn-1):
            lis_a.append(0.5*(lis_h[i]+lis_h[i+1])*self.mat_dat[i][0])
            lis_p.append(((lis_h[i+1]-lis_h[i])**2+self.mat_dat[i][0]**2)**0.5)
            lis_R.append(lis_a[i]/lis_p[i])
        P=sum(lis_p)+self.mat_dat[nn-2][1]+self.h0
        R=sum(lis_R)
        suman=0
        for i in range(len(self.mat_dat)):
            suman=suman+(lis_p[i]*lis_R[i]**(1/3))*(self.mat_dat[i][2])   
        np=(suman)/(P*R**(1/3))
        A=sum(lis_a)
        Q=(A**(5/3)*self.sf**(1/2))/(P**(2/3)*np)
        V=Q/A
        T=max(lis_x)
        g=9.807
        F=(V)/((g*((A)/(T)))**0.5)
        if (F<1):
            tex="Subcrítico"
        elif (F==1):
            tex="Crítico"
        elif (F>1):
            tex="Supercrítico"
        print("================Rugosidad ponderada Yen =================")
        print("Caudal                 (Q)[m3/s] :",Q)
        print("Área hidráulica          (A)[m2] :",A)
        print("Perímetro mojado          (P)[m] :",P)
        print("Rugosidad ponderada         (np) :",np)
        print("Velocidad               (V)[m/s] :",V)
        print("Espejo de agua            (T)[m] :",T)
        print("Número de Froude             (F) :",F)
        print("Tipo de flujo                    :",tex)
        print("=========================================================")
        return [Q,A,P,np,V,T,F,lis_x,lis_y]
       

    #def plot(lis_x,lis_y,lis_col,tex):
    def show_section(self,lis_col=["cyan","red","magenta"],alpha=0.9):
        res=self.results(show_results="FALSE")
        lis_x=res[7]
        lis_y=res[8]
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        xx=max(lis_x)
        plt.xlim(-0.05*xx,xx+0.05*xx)
        limy=min(lis_y)
        plt.ylim(+limy+0.1*limy,0-0.1*limy)
        datxy=[]
        for i in range(len(lis_x)):
            datxy.append([lis_x[i],lis_y[i]])
        shape=Polygon(datxy,color=lis_col[0],alpha=alpha)
        plt.plot(lis_x,lis_y,"-",color=lis_col[1])
        plt.plot(lis_x,lis_y,"o",color=lis_col[2])
        ax1.add_patch(shape)
        plt.grid()
        plt.title("Rugosidad ponderada Yen")
        plt.show()

#mat_dat=[[1.5,1.4,0.014],
#         [3,0.5,0.015]]
#sf=0.0005 #Pendiente de la linea de energia
#h0=1.4 #Profundidad en el margen (h0)[m]
#n0=0.016 #Coeficiente de Manning en el margen (n0)
#nf=0.016
#rr=RP_Yen(mat_dat,sf,h0,n0,nf)
#rr.results()
#rr.show_section(lis_col=["blue","red","yellow"],alpha=0.1)

class DVelocidad_maxima:
    def __init__(self,Q,n,Vmax,sf,z):
        self.Q=Q
        self.n=n
        self.Vmax=Vmax
        self.sf=sf
        self.z=z
    def results(self,show_iteration="FALSE",show_results="TRUE"):
        A=self.Q/self.Vmax
        R=((self.Q*self.n)/(self.sf**(1/2)*A))**(3/2)
        P=(A)/(R)
        if (show_results=="TRUE"):
            print("Área hidráulica         (A)[m2] :",A)    
            print("Radio hidráulico         (R)[m] :",R)
            print("Perímetro mojado         (P)[m] :",P)
        #dis=P**2-4*(self.z-2*(self.z**2+1)**0.5)*(-A)
        #print(dis)
        def f(y):
            res=(self.z-2*(self.z**2+1)**0.5)*y**2+P*y-A
            return res
        def df(y):
            res=2*(self.z-2*(self.z**2+1)**0.5)*y+P
            return res
        tol=0.0001
        y0=0.5
        de=6
        #y=newton(f,df,y0,tol,de)
        ##############################
        error=500
        i=0
        if (show_iteration=="TRUE"):
            print("N0".center(5),"y0".center(12),"f(y0)".center(12),"f'(y0)".center(12),
                  "y1".center(12),"error".center(12))
            print("=========================================================================")
        while (error>tol):
            i=i+1
            if (i>20):
                break
            y1=y0-(f(y0))/(df(y0))
            error=abs(y1-y0)
            if (show_iteration=="TRUE"):
                print(str(i).center(5),str(round(y0,de)).center(12),
                      str(round(f(y0),de)).center(12),
                      str(round(df(y0),de)).center(12),
                      str(round(y1,de)).center(12),str(round(error,de)).center(12))     
            y0=y1
        if (show_iteration=="TRUE"):
            print("=========================================================================")
        y=y1
        #######################################
        bb=P-2*y*(self.z**2+1)**0.5
        if (show_results=="TRUE"):
            print("Tirante normal          (yn)[m] :",y)
            print("Ancho de solera          (b)[m] :",bb)
        return [y,bb]
    def show_section(self,outline_color="silver",fill_color="cyan",width_size=7,height_size=5):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        b=res[1]
        z1=self.z
        z2=self.z
        T=z1*y+b+z2*y
        fig1=plt.figure(figsize=(width_size,height_size)) 
        ax1=fig1.add_subplot(xticks=[], yticks=[])
        dx1=-b/2-1.25*y*z1 
        dx2=b/2+1.25*y*z2
        plt.xlim(dx1-0.25*T,dx2+0.25*T)
        plt.ylim(-0.7*T+0.5*y,0.7*T+0.5*y)
        yyy=0
        xxx=0
        array=np.array([[-b/2+xxx,0+yyy],[-b/2-z1*y+xxx,y+yyy],
                          [b/2+z2*y+xxx,y+yyy],[b/2+xxx,0+yyy]])
        shape=Polygon(array,color=fill_color)
        ax1.add_patch(shape)
        plt.plot([-b/2-1.25*y*z1+xxx,-b/2+xxx,b/2+xxx,b/2+1.25*y*z2+xxx],
                   [1.25*y+yyy,0+yyy,0+yyy,1.25*y+yyy],"-",color=outline_color,lw=6)
        dy=0.05*abs(-dx1+dx2)
        if (b!=0):
            plt.plot([-b/2,b/2],[-dy,-dy],"-",alpha=1,color="black",lw=1)
            plt.plot([-b/2,-b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.plot([b/2,b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.text(0,-3*dy,str(round(b,3))+"  m",style="italic")
            plt.text(0-0.25*b,-3*dy,"b=",style="italic")

        plt.plot([b/2+z2*y+dy,b/2+z2*y+dy],[0,y],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y+dy-0.3*dy ,b/2+z2*y+dy+0.3*dy],[y,y],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y+dy-0.3*dy ,b/2+z2*y+dy+0.3*dy],[0,0],"-",alpha=1,color="black",lw=1)
        plt.text(b/2+z2*y+dy+0.3*dy,0.5*y,"y="+str(round(y,3))+"  m",style="italic")
            
        plt.plot([-b/2-z1*y,b/2+z2*y],[y+dy,y+dy],"-",alpha=1,color="black",lw=1)
        plt.plot([-b/2-z1*y,-b/2-z1*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y,b/2+z2*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.text(0,y+2*dy,str(round(T,3))+"  m",style="italic")
        plt.text(0-0.1*T,y+2*dy,"T=",style="italic")

        #if (b==0):
        plt.title("Canal no revestido, Velocidad máxima",style="italic",fontweight="bold")
        #elif (z1==0 and z2==0):
        #    plt.title("Tirante normal sección rectangular",style="italic",fontweight="bold")
        #elif (z1!=0 or z2!=0):
        #    plt.title("Tirante normal sección trapezoidal",style="italic",fontweight="bold") 
        plt.show()
    def show_section3d(self,lis_col=["red","lightcoral","cyan","blue"]):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        y=res[0]
        b=res[1]
        z1=self.z
        z2=self.z
        L=10
        rr=0.25
        h=y+rr*y
        points = np.array([[-0.5*L,-0.5*b,0],
                          [0.5*L, -0.5*b, 0 ],
                          [0.5*L,-0.5*b-z1*h, h],
                          [-0.5*L,-0.5*b-z1*h,h],
                          [-0.5*L,0.5*b,0],
                          [0.5*L, 0.5*b, 0 ],
                          [0.5*L,0.5*b+z2*h, h],
                          [-0.5*L,0.5*b+z2*h,h] 
                           ])
        #ancho=(z1*h+b+z2*h)/15
        ancho=0.2*(z1*y+b+z2*y)
        anchov=0.2*(y)
        points2 = np.array([[-0.5*L,-0.5*b-ancho-z1*ancho,-anchov],
                          [0.5*L, -0.5*b-ancho-z1*ancho, -anchov],
                          [0.5*L,-0.5*b-z1*h-ancho, h],
                          [-0.5*L,-0.5*b-z1*h-ancho,h],
                          [-0.5*L,0.5*b+ancho+z2*ancho,-anchov],
                          [0.5*L, 0.5*b+ancho+z2*ancho, -anchov],
                          [0.5*L,0.5*b+z2*h+ancho, h],
                          [-0.5*L,0.5*b+z2*h+ancho,h] 
                           ])
    
        points1 = np.array([[-0.5*L,-0.5*b,0],
                          [0.5*L, -0.5*b, 0 ],
                          [0.5*L,-0.5*b-z1*y, y],
                          [-0.5*L,-0.5*b-z1*y,y],
                          [-0.5*L,0.5*b,0],
                          [0.5*L, 0.5*b, 0 ],
                          [0.5*L,0.5*b+z2*y, y],
                          [-0.5*L,0.5*b+z2*y,y] 
                           ])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        r = [-5,5]
        X, Y = np.meshgrid(r,r)
        ax.set(xlim3d=(-0.5*L,0.5*L), xlabel='L')
        ax.set(ylim3d=(-0.5*b-z1*h-ancho-0.1*h,0.5*b+z2*h+ancho+0.1*h ),
               ylabel='')
        ax.set(zlim3d=(-ancho-0.3*h,h+0.3*h), zlabel='y')
        # vértices 
        #ax.scatter3D(points[:,0],points[:,1],points[:,2]) 
        #Lista de polígonos de los lados de la figura
        verts=[[points[0],points[1],points[2],points[3]],
               [points[4],points[5],points[6],points[7]],
               [points[0],points[1],points[5],points[4]]
               ]
        verts2=[[points2[0],points2[1],points2[2],points2[3]],
               [points2[4],points2[5],points2[6],points2[7]],
               [points2[0],points2[1],points2[5],points2[4]],
               [points[2],points[3],points2[3],points2[2]],
               [points[6],points[7],points2[7],points2[6]]  
               ]
        verts1=[[points1[0],points1[1],points1[2],points1[3]],
                [points1[4],points1[5],points1[6],points1[7]],
                [points1[0],points1[1],points1[5],points1[4]],
                [points1[0],points1[4],points1[7],points1[3]],
                [points1[1],points1[5],points1[6],points1[2]],
                [points1[2],points1[6],points1[7],points1[3]]
                ]
    #lados de la trama
        ax.add_collection3d(Poly3DCollection(verts,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))

        ax.add_collection3d(Poly3DCollection(verts2,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))
        ax.add_collection3d(Poly3DCollection(verts1,
                                             facecolors=lis_col[2],linewidths=1,
                                             edgecolors=lis_col[3], alpha=0.5))
        #if (self.b==0):
        #    plt.title("Tirante normal sección triangular")
        #elif (self.z1==0 and self.z2==0):
        plt.title("Canal no revestido, Velocidad máxima")
        #elif (self.z1!=0 or self.z2!=0):
        #    plt.title("Tirante normal sección trapezoidal")    
        plt.show()


#Q=60 #Caudal Q[m3/s]
#n=0.0222225348104 #Coeficiente de Manning
#Vmax=1.7 #Velocidad maxima
#sf=0.001 #Pendiente de la linea de energia
#z=2
#rr=DVelocidad_maxima(Q,n,Vmax,sf,z)
#rr.results(show_iteration="TRUE",show_results="TRUE")
#rr.show_section3d()

class DFuerza_tractiva:
    def __init__(self,Q,sf,n,d75,z,by,alfa,gamma,taup):
        self.Q=Q
        self.sf=sf
        self.n=n
        self.d75=d75
        self.z=z
        self.by=by
        self.alfa=alfa
        self.gamma=gamma
        self.taup=taup
    def results(self,show_results="TRUE"):
        Gamma=atan(1/self.z)
        K=(1-(sin(Gamma))**2/(sin((pi/180)*self.alfa))**2)**0.5
        #Calculo del esfuerzo tangencial critico en los taludes
        es=epsilons(self.z,self.by)
        ep=epsilonp(self.z,self.by)
        taus=K*self.taup
        ys=(taus)/(es*self.gamma*self.sf)
        yp=(self.taup)/(ep*self.gamma*self.sf)
        y=min([ys,yp])
        b=self.by*y
        A=self.z*y**2+b*y
        P=2*y*(self.z**2+1)**0.5+b
        q=(A**(5/3)*self.sf**(1/2))/(P**(2/3)*self.n)
        V=(q)/(A)
        g=9.807
        T=2*self.z*y+b
        F=(V)/((g*((A)/(T)))**0.5)
        ########
        if (show_results=="TRUE"):
            print("Angulo del los lados del canal                   (Г)[°] :",(180/pi)*Gamma)
            print("Angulo de reposo                                 (α)[°] :",self.alfa)
            print("Angulo del los lados del canal                 (Г)[rad] :",Gamma)
            print("Angulo de reposo                               (α)[rad] :",(pi/180)*self.alfa)
            print("Razon de fuerzas tractivas                          (k) :",K)
            print("(Factor de corrección en el fondo                  (εp) :",ep)
            print("(Factor de corrección en el talud                  (εs) :",es)
            print("Esfuerzo de corte permisible en el fondo    (τp)[kg/m2] :",self.taup)
            print("Esfuerzo de corte permisible en los taludes (τs)[kg/m2] :",taus)
            print("Tirante                                         (yp)[m] :",yp)
            print("Tirante                                         (ys)[m] :",ys)
            print("Ancho de solera                                  (b)[m] :",b)
            print("Area hidráulica                                 (A)[m2] :",A)
            print("Perímetro mojado                                 (P)[m] :",P)
            print("Espejo de agua                                   (T)[m] :",T)
            print("Velocidad                                      (v)[m/s] :",V)
            print("Número de Froude                                    (F) :",F)
            print("Caudal de diseño                             (Qd)[m3/s] :",self.Q)
            print("Caudal calculado                             (Qc)[m3/s] :",q)
        return [y,b]
    def show_section(self,outline_color="silver",fill_color="cyan",width_size=7,height_size=5):
        res=self.results(show_results="FALSE")
        y=res[0]
        b=res[1]
        z1=self.z
        z2=self.z
        T=z1*y+b+z2*y
        fig1=plt.figure(figsize=(width_size,height_size)) 
        ax1=fig1.add_subplot(xticks=[], yticks=[])
        dx1=-b/2-1.25*y*z1 
        dx2=b/2+1.25*y*z2
        plt.xlim(dx1-0.25*T,dx2+0.25*T)
        plt.ylim(-0.7*T+0.5*y,0.7*T+0.5*y)
        yyy=0
        xxx=0
        array=np.array([[-b/2+xxx,0+yyy],[-b/2-z1*y+xxx,y+yyy],
                          [b/2+z2*y+xxx,y+yyy],[b/2+xxx,0+yyy]])
        shape=Polygon(array,color=fill_color)
        ax1.add_patch(shape)
        plt.plot([-b/2-1.25*y*z1+xxx,-b/2+xxx,b/2+xxx,b/2+1.25*y*z2+xxx],
                   [1.25*y+yyy,0+yyy,0+yyy,1.25*y+yyy],"-",color=outline_color,lw=6)
        dy=0.05*abs(-dx1+dx2)
        if (b!=0):
            plt.plot([-b/2,b/2],[-dy,-dy],"-",alpha=1,color="black",lw=1)
            plt.plot([-b/2,-b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.plot([b/2,b/2],[-dy-0.3*dy,-dy+0.3*dy],"-",alpha=1,color="black",lw=1)
            plt.text(0,-3*dy,str(round(b,3))+"  m",style="italic")
            plt.text(0-0.25*b,-3*dy,"b=",style="italic")

        plt.plot([b/2+z2*y+dy,b/2+z2*y+dy],[0,y],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y+dy-0.3*dy ,b/2+z2*y+dy+0.3*dy],[y,y],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y+dy-0.3*dy ,b/2+z2*y+dy+0.3*dy],[0,0],"-",alpha=1,color="black",lw=1)
        plt.text(b/2+z2*y+dy+0.3*dy,0.5*y,"y="+str(round(y,3))+"  m",style="italic")
            
        plt.plot([-b/2-z1*y,b/2+z2*y],[y+dy,y+dy],"-",alpha=1,color="black",lw=1)
        plt.plot([-b/2-z1*y,-b/2-z1*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.plot([b/2+z2*y,b/2+z2*y],[y+dy-0.3*dy,y+dy+0.3*dy],"-",alpha=1,color="black",lw=1)
        plt.text(0,y+2*dy,str(round(T,3))+"  m",style="italic")
        plt.text(0-0.1*T,y+2*dy,"T=",style="italic")

        #if (b==0):
        plt.title("Canal no revestido, Fuerza tractiva",style="italic",fontweight="bold")
        #elif (z1==0 and z2==0):
        #    plt.title("Tirante normal sección rectangular",style="italic",fontweight="bold")
        #elif (z1!=0 or z2!=0):
        #    plt.title("Tirante normal sección trapezoidal",style="italic",fontweight="bold") 
        plt.show()
    def show_section3d(self,lis_col=["red","lightcoral","cyan","blue"]):
        res=self.results(show_results="FALSE")
        y=res[0]
        b=res[1]
        z1=self.z
        z2=self.z
        L=10
        rr=0.25
        h=y+rr*y
        points = np.array([[-0.5*L,-0.5*b,0],
                          [0.5*L, -0.5*b, 0 ],
                          [0.5*L,-0.5*b-z1*h, h],
                          [-0.5*L,-0.5*b-z1*h,h],
                          [-0.5*L,0.5*b,0],
                          [0.5*L, 0.5*b, 0 ],
                          [0.5*L,0.5*b+z2*h, h],
                          [-0.5*L,0.5*b+z2*h,h] 
                           ])
        #ancho=(z1*h+b+z2*h)/15
        ancho=0.2*(z1*y+b+z2*y)
        anchov=0.2*(y)
        points2 = np.array([[-0.5*L,-0.5*b-ancho-z1*ancho,-anchov],
                          [0.5*L, -0.5*b-ancho-z1*ancho, -anchov],
                          [0.5*L,-0.5*b-z1*h-ancho, h],
                          [-0.5*L,-0.5*b-z1*h-ancho,h],
                          [-0.5*L,0.5*b+ancho+z2*ancho,-anchov],
                          [0.5*L, 0.5*b+ancho+z2*ancho, -anchov],
                          [0.5*L,0.5*b+z2*h+ancho, h],
                          [-0.5*L,0.5*b+z2*h+ancho,h] 
                           ])
    
        points1 = np.array([[-0.5*L,-0.5*b,0],
                          [0.5*L, -0.5*b, 0 ],
                          [0.5*L,-0.5*b-z1*y, y],
                          [-0.5*L,-0.5*b-z1*y,y],
                          [-0.5*L,0.5*b,0],
                          [0.5*L, 0.5*b, 0 ],
                          [0.5*L,0.5*b+z2*y, y],
                          [-0.5*L,0.5*b+z2*y,y] 
                           ])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        r = [-5,5]
        X, Y = np.meshgrid(r,r)
        ax.set(xlim3d=(-0.5*L,0.5*L), xlabel='L')
        ax.set(ylim3d=(-0.5*b-z1*h-ancho-0.1*h,0.5*b+z2*h+ancho+0.1*h ),
               ylabel='')
        ax.set(zlim3d=(-ancho-0.3*h,h+0.3*h), zlabel='y')
        # vértices 
        #ax.scatter3D(points[:,0],points[:,1],points[:,2]) 
        #Lista de polígonos de los lados de la figura
        verts=[[points[0],points[1],points[2],points[3]],
               [points[4],points[5],points[6],points[7]],
               [points[0],points[1],points[5],points[4]]
               ]
        verts2=[[points2[0],points2[1],points2[2],points2[3]],
               [points2[4],points2[5],points2[6],points2[7]],
               [points2[0],points2[1],points2[5],points2[4]],
               [points[2],points[3],points2[3],points2[2]],
               [points[6],points[7],points2[7],points2[6]]  
               ]
        verts1=[[points1[0],points1[1],points1[2],points1[3]],
                [points1[4],points1[5],points1[6],points1[7]],
                [points1[0],points1[1],points1[5],points1[4]],
                [points1[0],points1[4],points1[7],points1[3]],
                [points1[1],points1[5],points1[6],points1[2]],
                [points1[2],points1[6],points1[7],points1[3]]
                ]
    #lados de la trama
        ax.add_collection3d(Poly3DCollection(verts,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))

        ax.add_collection3d(Poly3DCollection(verts2,
                                             facecolors=lis_col[0],linewidths=1,
                                             edgecolors=lis_col[1], alpha=0.1))
        ax.add_collection3d(Poly3DCollection(verts1,
                                             facecolors=lis_col[2],linewidths=1,
                                             edgecolors=lis_col[3], alpha=0.5))
        #if (self.b==0):
        #    plt.title("Tirante normal sección triangular")
        #elif (self.z1==0 and self.z2==0):
        plt.title("Canal no revestido, Fuerza tractiva")
        #elif (self.z1!=0 or self.z2!=0):
        #    plt.title("Tirante normal sección trapezoidal")    
        plt.show()
        
#Q=6
#sf=0.001
#n=0.03
#d75=40
#z=1.5#talud(z=0,1.5,2)
#Cs=0.90#coeficiente de sinousidad
#by=0.05#relacion b entre y
#alfa=37.5#angulo de reposo [grados sexagesimales]
#gamma=1000#[kg/m3]
#taup=(d75)/(13)#kg/m2
#################################################
#rr=DFuerza_tractiva(Q,sf,n,d75,z,by,alfa,gamma,taup)
#rr.results(show_results="TRUE")
#rr.show_section3d()





