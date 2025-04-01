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
def EET(Q,b,z1,z2,alfa,y):
    g=9.807
    A=0.5*(z1+z2)*y**2+b*y
    V=Q/A
    res=y+alfa*(V**2)/(2*g)
    return res
def show_EYT(b,z1,z2,Q,alfa,lis_py=[],lis_col=["blue","red","orange","gray"],posx=0.80,posy=0.5,tam=0.09,dimy=3*0.6):
    nn=50
    ##########################
    error=500
    tol=0.001
    g=9.807
    yc0=0.5
    yc1=1
    while (error>tol):
        A=b*yc0+0.5*(z1+z2)*yc0**2
        T=(z1+z2)*yc0+b
        dA=b+(z1+z2)*yc0
        dT=z1+z2
        fy=(A**(3))/(T)-(Q**2)/(g/alfa)
        dfy=(3*A**2*T*dA-A**3*dT)/(T**2)
        yc1=yc0-fy/dfy
        error=abs(yc1-yc0)
        yc0=yc1
    ##################
    yc=yc1
    fac=tam
    inc=(dimy)/nn
    Y=[inc*(i+1) for i in range(nn)]
    E=[EET(Q,b,z1,z2,alfa,Y[i]) for i in range(len(Y))]    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, aspect='equal')
    plt.xlim(0,1.5*max(Y))
    plt.ylim(0,1.1*max(Y))
    yyy=posy*(dimy)
    xxx=posx*(1.5*max(Y))
    fac1=(dimy)/(yc)
    yyc=fac*yc*fac1
    bb=fac*b*fac1
    array=np.array([[-bb/2+xxx,0+yyy],[-bb/2-z1*yyc+xxx,yyc+yyy],
                    [bb/2+z2*yyc+xxx,yyc+yyy],[bb/2+xxx,0+yyy]])
    shape=Polygon(array,color=lis_col[2])
    ax1.add_patch(shape)
    plt.plot([-bb/2-1.25*yyc*z1+xxx,-bb/2+xxx,bb/2+xxx,bb/2+1.25*yyc*z2+xxx],
            [1.25*yyc+yyy,0+yyy,0+yyy,1.25*yyc+yyy],"-",color=lis_col[3])                    
    plt.plot(Y,Y,"-",color="silver",lw=2)
    plt.plot(E,Y,"-",color=lis_col[0],lw=2,label="Y(E)")
    for i in range(len(lis_py)):
        lista=[lis_py[i],EET(Q,b,z1,z2,alfa,lis_py[i])]
        plt.plot([0,lista[1]],[lista[0],lista[0]],"--",color=lis_col[1])
        plt.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_col[1])
        plt.plot(lista[1],lista[0],"mo",markersize=11,mec="k",alpha=1,mew=1.5)
    plt. xlabel ("Energía específica (E)",style="italic",fontweight="bold")
    plt. ylabel ("Tirante (Y)",style="italic",fontweight="bold")
    plt.title("E vs Y",style="italic",fontweight="bold")
    #nx = 1
    plt.legend(loc="upper right")
    plt.grid()
    plt.show()













class Tirantec_t:
    def __init__(self,Q,b,z1,z2,alfa):
        self.Q=Q
        self.b=b
        self.z1=z1
        self.z2=z2
        self.alfa=alfa
    def results(self,yc0=0.5,decimal_number=6,tol=0.0001,show_iteration="FALSE",show_results="TRUE"):
        error=50
        if (show_iteration=="TRUE"):
            print("yi[m]".center(10),"A[m2]".center(10),
                  "dA".center(10),"T".center(10),"dT".center(10),
                  "f(y)".center(10),"df(y)".center(10),"yi+1[m]".center(10),
                  "error".center(10))
            print(" ==================================================================================================")
        while (error>tol):
            g=9.807
            A=self.b*yc0+0.5*(self.z1+self.z2)*yc0**2
            T=(self.z1+self.z2)*yc0+self.b
            dA=self.b+(self.z1+self.z2)*yc0
            dT=self.z1+self.z2
            fy=(A**(3))/(T)-(self.Q**2)/(g/self.alfa)
            dfy=(3*A**2*T*dA-A**3*dT)/(T**2)
            yc1=yc0-fy/dfy
            error=abs(yc1-yc0)
            if (show_iteration=="TRUE"):
                print(str(round(yc0,6)).center(10),str(round(A,6)).center(10),
                      str(round(dA,6)).center(10),str(round(T,6)).center(10),str(round(dT,6)).center(10),
                      str(round(fy,6)).center(10),str(round(dfy,6)).center(10),str(round(yc1,6)).center(10),
                      str(round(error,6)).center(10))
            yc0=yc1
            if (show_iteration=="TRUE"):
                print(" ==================================================================================================")    
        g=9.807
        P=((1+self.z1**2)**0.5+(1+self.z2**2)**0.5)*yc1+self.b
        R=A/P
        V=self.Q/A
        F=(V)/((g*(A)/(T))**0.5)
        E=yc1+(V**2)/(2*g)
        if (show_results=="TRUE"):
            print(" Tirante crítico          yc[m] :",yc1)
            print(" Área hidráulica         A[m2] :",A)
            print(" Espejo de agua           T[m] :",T)
            print(" Número de Froude            F :",1)
            print(" Perímetro mojado         P[m] :",P)
            print(" Radio hidráulico         R[m] :",R)
            print(" Velocidad              V[m/s] :",V)
            print(" Energía específica E[m-kg/kg] :",E)
        return [yc1,A,T,1,P,R,V,E]
        

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
            plt.title("Tirante crítico sección triangular")
        elif (self.z1==0 and self.z2==0):
            plt.title("Tirante crítico sección rectangular")
        elif (self.z1!=0 or self.z2!=0):
            plt.title("Tirante crítico sección trapezoidal")    
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
            plt.title("Tirante crítico sección triangular",style="italic",fontweight="bold")
        elif (self.z1==0 and self.z2==0):
            plt.title("Tirante crítico sección rectangular",style="italic",fontweight="bold")
        elif (self.z1!=0 or self.z2!=0):
            plt.title("Tirante crítico sección trapezoidal",style="italic",fontweight="bold") 
        plt.show()
    def show_EY(self,lis_py=[],lis_col=["blue","red","orange","gray"],posx=0.80,posy=0.5,tam=0.09,dimy=3*0.6):
        nn=50
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        yc=res[0]
        fac=tam
        inc=(dimy)/nn
        Y=[inc*(i+1) for i in range(nn)]
        E=[EET(self.Q,self.b,self.z1,self.z2,self.alfa,Y[i]) for i in range(len(Y))]    
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.xlim(0,1.5*max(Y))
        plt.ylim(0,1.1*max(Y))
        yyy=posy*(dimy)
        xxx=posx*(1.5*max(Y))
        fac1=(dimy)/(yc)
        yc=fac*yc*fac1
        b=fac*self.b*fac1
        array=np.array([[-b/2+xxx,0+yyy],[-b/2-self.z1*yc+xxx,yc+yyy],
                        [b/2+self.z2*yc+xxx,yc+yyy],[b/2+xxx,0+yyy]])
        shape=Polygon(array,color=lis_col[2])
        ax1.add_patch(shape)
        plt.plot([-b/2-1.25*yc*self.z1+xxx,-b/2+xxx,b/2+xxx,b/2+1.25*yc*self.z2+xxx],
                 [1.25*yc+yyy,0+yyy,0+yyy,1.25*yc+yyy],"-",color=lis_col[3])                    
        plt.plot(Y,Y,"-",color="silver",lw=2)
        plt.plot(E,Y,"-",color=lis_col[0],lw=2,label="Y(E)")
        for i in range(len(lis_py)):
            lista=[lis_py[i],EET(self.Q,self.b,self.z1,self.z2,self.alfa,lis_py[i])]
            plt.plot([0,lista[1]],[lista[0],lista[0]],"--",color=lis_col[1])
            plt.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_col[1])
            plt.plot(lista[1],lista[0],"mo",markersize=11,mec="k",alpha=1,mew=1.5)
        plt. xlabel ("Energía específica (E)",style="italic",fontweight="bold")
        plt. ylabel ("Tirante (Y)",style="italic",fontweight="bold")
        plt.title("E vs Y",style="italic",fontweight="bold")
        #nx = 1
        plt.legend(loc="center right")
        plt.grid()
        plt.show()
    def show_EYQ(self,lis_colQ=["yellow","orange","red","gray","blue","cyan","gray"],lis_Q=[0.5,0.6],dimy=1*0.6):
        nn=50
        inc=(dimy)/nn
        Y=[inc*(i+1) for i in range(nn)]
        lis_E=[]
        for i in range(len(lis_Q)):
            E=[EET(lis_Q[i],self.b,self.z1,self.z2,self.alfa,Y[j]) for j in range(nn)]
            lis_E.append(E)    
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.xlim(0,1.5*max(Y))
        plt.ylim(0,1.1*max(Y))
        plt.plot(Y,Y,"-",color="silver",lw=2)
        for i in range(len(lis_Q)):
            plt.plot(lis_E[i],Y,"-",color=lis_colQ[i],lw=2,label="Q[m3/s]="+str(round(lis_Q[i],3)))
        plt. xlabel ("Energía específica (E)",style="italic",fontweight="bold")
        plt. ylabel ("Tirante (Y)",style="italic",fontweight="bold")
        plt.title("E vs Y",style="italic",fontweight="bold")
        #nx = 1
        plt.legend(loc="upper left")
        plt.grid()
        plt.show()

#Q=5      #Caudal Q[m3/s]
#b=1      #Ancho de solera b[m]
#z1=2     #Talud 1
#z2=2     #Talud 2
#alfa=1
#rr=Tirantec_t(Q,b,z1,z2,alfa)
#rrr=rr.results(yc0=0.5,decimal_number=6,tol=0.0001,show_iteration="FALSE",show_results="TRUE")
#yc=rrr[0]
#rr.show_section()
#rr.show_EY(lis_py=[yc,0.1,0.5,0.2,0.6],lis_col=["red","blue","cyan","yellow"],dimy=2*yc,posx=0.5,posy=0.5,tam=0.5)
#rr.show_EYQ(lis_Q=[0.5,0.6,0.7,0.8,0.9],lis_colQ=["blue","red","orange","yellow","cyan"])


def EEC(Q,D,y):
    g=9.807
    A=0.25*D**2*acos(1-2*y/D)-0.5*D**2*(y/D-y**2/D**2)**0.5*(1-2*y/D)
    V=Q/A
    res=y+(V**2)/(2*g)
    return res
class Tirantec_c:
    def __init__(self,Q,D,alfa):
        self.Q=Q
        self.D=D
        self.alfa=alfa
    def results(self,decimal_number=6,tol=0.0001,show_iteration="FALSE",show_results="TRUE"):
        g=9.807
        error=50
        if (show_iteration=="TRUE"):
            print("yi[m]".center(10),"A[m2]".center(10),
                  "dA".center(10),"T".center(10),"dT".center(10),
                  "f(y)".center(10),"df(y)".center(10),"yi+1[m]".center(10),
                  "error".center(10))
            print(" ==================================================================================================")
        lis_Cau=[]
        for i in range(4):
            yc0=(i+1)*self.D/4
            a=0.25*self.D**2*acos(1-2*yc0/self.D)-0.5*self.D**2*(yc0/self.D-yc0**2/self.D**2)**0.5*(1-2*yc0/self.D)
            p=self.D*acos(1-2*yc0/self.D)
            r=a/p
            teta=2*acos(1-(2*yc0)/self.D)
            tt=self.D*sin(teta/2)
            Cau=((a**3*g)/(tt))**0.5
            lis_Cau.append(Cau)
        if (0<self.Q and  self.Q <=lis_Cau[0]):
            yc0=0*self.D/4+self.D/8
        elif (lis_Cau[0]<self.Q and self.Q<=lis_Cau[1]):
            yc0=1*self.D/4+self.D/8
        elif (lis_Cau[1]<self.Q):
            yc0=0.75*self.D    
        while (error>tol):
            A=0.25*self.D**2*acos(1-2*yc0/self.D)-0.5*self.D**2*(yc0/self.D-yc0**2/self.D**2)**0.5*(1-2*yc0/self.D)
            T=self.D*sin(acos(1-(2*yc0)/self.D))
            dA=0.25*self.D**2*((2/self.D)/(1-(1-2*yc0/self.D)**2)**0.5)-0.5*self.D**2*((1-2*yc0/self.D)*((1/self.D-2*yc0/self.D**2)/(2*(yc0/self.D-yc0**2/self.D**2)**0.5))+(yc0/self.D-yc0**2/self.D**2)**0.5*(-2/self.D))
            dT=(2*cos(acos(1-(2*yc0)/self.D)))*(1/(1-(1-(2*yc0)/self.D)**2)**0.5)
            fy=(A**(3))/(T)-(self.Q**2)/(g/self.alfa)
            dfy=(3*A**2*T*dA-A**3*dT)/(T**2)
            yc1=yc0-fy/dfy
            error=abs(yc1-yc0)
            if (show_iteration=="TRUE"):
                print(str(round(yc0,decimal_number)).center(10),str(round(A,decimal_number)).center(10),
                      str(round(dA,decimal_number)).center(10),str(round(T,decimal_number)).center(10),str(round(dT,decimal_number)).center(10),
                      str(round(fy,decimal_number)).center(10),str(round(dfy,decimal_number)).center(10),str(round(yc1,decimal_number)).center(10),
                      str(round(error,decimal_number)).center(10))
            yc0=yc1
        if (show_iteration=="TRUE"):
            print(" ==================================================================================================")    
        g=9.807
        P=self.D*acos(1-2*yc1/self.D)
        R=A/P
        V=self.Q/A
        F=(V)/((g*(A)/(T))**0.5)
        E=yc1+self.alfa*(V**2)/(2*g)
        if (show_results=="TRUE"):
            print(" Tirante crítico         yc[m] :",yc1)
            print(" Área hidráulica         A[m2] :",A)
            print(" Espejo de agua           T[m] :",T)
            print(" Número de Froude            F :",1)
            print(" Perímetro mojado         P[m] :",P)
            print(" Radio hidráulico         R[m] :",R)
            print(" Velocidad              V[m/s] :",V)
            print(" Energía específica E[m-kg/kg] :",E)
        return [yc1,A,T,1,P,R,V,E]
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
        plt.title("Tirante crítico sección circular",style="italic",fontweight="bold")
    
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
        #lis_QY=lisQY(self.D,self.n,self.sf,N)
        #lis_VY=lisVY(self.D,self.n,self.sf,N)
        #QQ=lis_QY[0]
        #VV=lis_VY[0]     
        #Y=lis_QY[1]
        #Qmax=max(QQ)
        #Vmax=max(VV)
        #plt.text(-self.D/2-self.D/2.5,-self.D/2-self.D/5,"Qmax="+str(round(Qmax,3))+" m3/s",style="italic")
        #plt.text(-self.D/2-self.D/2.5,-self.D/2-self.D/3.5,"Vmax="+str(round(Vmax,3))+" m/s",style="italic")
        plt.show()
    def show_EY(self,lis_py=[],lis_col=["blue","red","orange","gray"],posx=0.80,posy=0.5,tam=0.4):
        res=self.results(show_iteration="FALSE",show_results="FALSE")
        nn=50
        yc=res[0]
        fac=tam
        inc=(self.D)/nn
        Y=[inc*(i+1) for i in range(nn)]
        E=[EEC(self.Q,self.D,Y[i]) for i in range(len(Y))]        
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.xlim(0,1.5*max(Y))
        plt.ylim(0,1.1*max(Y))
        yyy=posy*self.D
        xxx=posx*(1.5*max(Y))
    #fac1=(1.1*max(Y))/(D)
        yc=fac*yc
        D=fac*self.D
#################
        y=yc
        thetar=2*acos(1-(2*y)/(D))
        thetag=(180/pi)*thetar
        theta11=180+(90-thetag/2)
        theta22=270+(thetag/2)
        ax1.add_patch(patches.Wedge(center=(0+xxx,0+yyy), r=D/2+D/35, theta1=0
                    , theta2=360, facecolor=lis_col[3], label="Test"))
        ax1.add_patch(patches.Wedge(center=(0+xxx,0+yyy), r=D/2, theta1=0
                    , theta2=360, facecolor="white", label="Test"))
        ax1.add_patch(patches.Wedge(center=(0+xxx,0+yyy), r=D/2, theta1=theta11
                    , theta2=theta22, facecolor=lis_col[2], label="Test"))
        array=np.array([[0+xxx,0+yyy],[-D/2*sin(thetar/2)+xxx,-(D/2-y)+yyy],[D/2*sin(thetar/2)+xxx,-(D/2-y)+yyy]])
        if (y<=D/2):
            shape=Polygon(array,color="white")
        elif (y>D/2):
            shape=Polygon(array,color=lis_col[2])
        ax1.add_patch(shape)    
#####################               
        plt.plot(Y,Y,"-",color="silver",lw=2)
        plt.plot(E,Y,"-",color=lis_col[0],lw=2,label="Y(E)")
        for i in range(len(lis_py)):
            lista=[lis_py[i],EEC(self.Q,self.D,lis_py[i])]
            plt.plot([0,lista[1]],[lista[0],lista[0]],"--",color=lis_col[1])
            plt.plot([lista[1],lista[1]],[0,lista[0]],"--",color=lis_col[1])
            plt.plot(lista[1],lista[0], "mo",markersize=11,mec="k",alpha=1,mew=1.5)
        plt. xlabel("Energía específica (E)",style="italic",fontweight="bold")
        plt. ylabel("Tirante (Y)",style="italic",fontweight="bold")
        plt.title("E vs Y",style="italic",fontweight="bold")
        #nx = 1
        plt.grid()
        #plt.legend(loc="upper left")
        plt.show()
    

#Q=5
#D=2.5
#tol=0.0001
######################
#alfa=1
#rr=Tirantec_c(Q,D,alfa)
#rrr=rr.results()
#rr.show_section()
#yc=rrr[0]
#rr.show_EY(lis_py=[yc])

def Ey1(Q,y1,H1,H2,z1,z2,b1,b2,y2,K):
    g=9.807
    A1=z1*y1**2+b1*y1
    A2=z2*y2**2+b2*y2
    res=-(H2-H1)+y2+(Q**2)/(2*g*A2**2)+K*((Q**2)/(2*g*A1**2)-(Q**2)/(2*g*A2**2))-y1-(Q**2)/(2*g*A1**2)
    return res
def Ey1e(Q,y1,H1,H2,z1,z2,b1,b2,y2,K):
    g=9.807
    A1=z1*y1**2+b1*y1
    A2=z2*y2**2+b2*y2
    res=(H1-H2)+y2+(Q**2)/(2*g*A2**2)+K*(-(Q**2)/(2*g*A1**2)+(Q**2)/(2*g*A2**2))-y1-(Q**2)/(2*g*A1**2)
    return res
def dEy1(Q,y1,H1,H2,z1,z2,b1,b2,y2,K):
    g=9.807
    A1=z1*y1**2+b1*y1
    dA1=2*z1*y1+b1
    res=-K*((Q**2*A1**(-3))/(g))*dA1-1+((Q**2*A1**(-3))/(g))*dA1
    return res

def dEy1e(Q,y1,H1,H2,z1,z2,b1,b2,y2,K):
    g=9.807
    A1=z1*y1**2+b1*y1
    dA1=2*z1*y1+b1
    res=K*((Q**2*A1**(-3))/(g))*dA1-1+((Q**2*A1**(-3))/(g))*dA1
    return res

def solvey1(Q,H1,H2,z1,z2,b1,b2,y2,K):
    error=500
    tol=0.0001
    y0=1
    while(error>tol):
        yf=y0-Ey1(Q,y0,H1,H2,z1,z2,b1,b2,y2,K)/dEy1(Q,y0,H1,H2,z1,z2,b1,b2,y2,K)
        error=abs(yf-y0)
        y0=yf
    return yf
def solvey1e(Q,H1,H2,z1,z2,b1,b2,y2,K):
    error=500
    tol=0.0001
    y0=1
    while(error>tol):
        yf=y0-Ey1e(Q,y0,H1,H2,z1,z2,b1,b2,y2,K)/dEy1e(Q,y0,H1,H2,z1,z2,b1,b2,y2,K)
        error=abs(yf-y0)
        y0=yf
    return yf 
def rev(lis):
    res=[]
    n=len(lis)
    for i in range(n):
        res.append(lis[n-1-i])
    return res
class Transicion_AS:
    def __init__(self,Q,bf,bc,zc,yc,HZ,K,nt):
        self.Q=Q
        self.bf=bf
        self.bc=bc
        self.zc=zc
        self.yc=yc
        self.HZ=HZ
        self.K=K
        self.nt=nt
    def results(self,decimal_number=6,show_results="TRUE"):
        zf=0
        b=0.5*(self.bc-self.bf)
        L=4.7*b+1.65*self.zc*self.yc
        L=int(L+1)
        nb=0.8-0.26*self.zc**0.5
        incx=L/self.nt
        x=[i*incx for i in range(self.nt+1)]
        lis_b=[]
        lis_z=[]
        for i in range(self.nt+1):
            lis_b.append(self.bf+(self.bc-self.bf)*(x[i]/L)*(1-(1-x[i]/L)**nb))
            lis_z.append(self.zc*(1-(1-x[i]/L)**0.5))
        lis_inch=[0]
        suma=0
        for i in range(self.nt):
            suma=suma+(self.HZ/L)*(x[i+1]-x[i])
            lis_inch.append(suma)
        lis_y1r=[self.yc]
        for i in range(self.nt):
            lis_y1r.append(solvey1(self.Q,lis_inch[self.nt-1-i],lis_inch[self.nt-i],
                                   lis_z[self.nt-1-i],lis_z[self.nt-i],
                                   lis_b[self.nt-1-i],lis_b[self.nt-i],lis_y1r[i],self.K))
        lis_y1=rev(lis_y1r)
        lis_v=[]
        lis_E=[]
        g=9.807
        for i in range(len(lis_b)):
            A=lis_z[i]*lis_y1[i]**2+lis_b[i]*lis_y1[i]
            lis_v.append(self.Q/A)
            lis_E.append(lis_y1[i]+(lis_v[i]**2)/(2*g))
        if (show_results=="TRUE"):    
            print("x[m]".center(12),"b[m]".center(12),"z".center(12),"HZ[m]".center(12),"y[m]".center(12),"v[m/s]".center(12),"E[m]".center(12))
            print("============================================================================================")

            for i in range(len(lis_v)):
                print(str(round(x[i],decimal_number)).center(12),str(round(lis_b[i],decimal_number)).center(12),
                      str(round(lis_z[i],decimal_number)).center(12),
                      str(round(lis_inch[i],decimal_number)).center(12),str(round(lis_y1[i],decimal_number)).center(12),
                      str(round(lis_v[i],decimal_number)).center(12),str(round(lis_E[i],decimal_number)).center(12))
        return [x,lis_b,lis_z,lis_inch,lis_y1,lis_v,lis_E,L,self.HZ]
    def show_transicion(self,lis_col=["cyan","black"]):
        lis_res=self.results(show_results="FALSE")
        lis_x=lis_res[0]
        lis_b=lis_res[1]
        lis_z=lis_res[2]
        lis_inch=lis_res[3]
        lis_y1=lis_res[4]
        lis_v=lis_res[5]
        lis_E=lis_res[6]
        L=lis_res[7]
        HZ=lis_res[8]
####################################
        bmin=min(lis_b)
        zmin=min(lis_z)
        ymin=min(lis_y1)
        bmax=max(lis_b)
        zmax=max(lis_z)
        ymax=max(lis_y1)
        desh=rev(lis_inch)
        points1 = np.array([
            [0,-0.5*bmin-zmin*ymin,ymin+HZ],
            [0,-0.5*bmin,0+HZ],
            [0,0.5*bmin,0+HZ],
            [0,0.5*bmin+zmin*ymin,ymin+HZ],
            [L,-0.5*bmax-zmax*ymax,ymax],
            [L,-0.5*bmax,0],
            [L,0.5*bmax, 0],
            [L,0.5*bmax+zmax*ymax,ymax]
            ])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        r = [-1,1]
        X, Y = np.meshgrid(r,r)
        # plot vertices
        ax.scatter3D(points1[:,0],points1[:,1],points1[:,2])
        verts1=[]
        for i in range(len(lis_b)):
            verts1.append([[lis_x[i],-0.5*lis_b[i]-lis_z[i]*lis_y1[i],lis_y1[i]+desh[i]],
                           [lis_x[i],-0.5*lis_b[i],0+desh[i]],
                           [lis_x[i],0.5*lis_b[i], 0+desh[i]],
                           [lis_x[i],0.5*lis_b[i]+lis_z[i]*lis_y1[i],lis_y1[i]+desh[i]]])
        ax.add_collection3d(Poly3DCollection(verts1,facecolors=lis_col[0],linewidths=1,edgecolors=lis_col[1], alpha=.25))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.title("Transición de salida (expansión)")
        plt.show()

#Q=5#Caudal [m3/s]
#bc=5#Ancho de solera seccion c [m]
#bf=3#Ancho de solera seccion f [m]
#zc=2#Talud en la seccion c
#yc=1.0207 
#HZ=0.1#Diferencia de cotas [m]
#K=0.2#Coeficiente de perdidas en la transicion
#nt=10
#rr=Transicion_AS(Q,bf,bc,zc,yc,HZ,K,nt)
#rr.results()
#rr.show_transicion()
class Transicion_AE:
    def __init__(self,Q,bf,bc,zf,yc,HZ,K,nt):
        self.Q=Q
        self.bf=bf
        self.bc=bc
        self.zf=zf
        self.yc=yc
        self.HZ=HZ
        self.K=K
        self.nt=nt
    def results(self,decimal_number=6,show_results="TRUE"):
        zc=0
        b=0.5*(self.bf-self.bc)
        #L=4.7*b+1.65*zf*yf
        L=100
        #L=int(L+1)
        nb=0.8-0.26*self.zf**0.5
        incx=L/self.nt
        x=[i*incx for i in range(self.nt+1)]
        lis_b=[]
        lis_z=[]
        for i in range(self.nt+1):
            lis_b.append(self.bc+(self.bf-self.bc)*(x[i]/L)*(1-(1-x[i]/L)**nb))
            lis_z.append(self.zf*(1-(1-x[i]/L)**0.5))
        lis_inch=[0]
        suma=0
        for i in range(self.nt):
            suma=suma+(self.HZ/L)*(x[i+1]-x[i])
            lis_inch.append(suma)
        lis_y1r=[self.yc]
        for i in range(self.nt):
            lis_y1r.append(solvey1e(self.Q,lis_inch[i+1],lis_inch[i],lis_z[i+1],lis_z[i],lis_b[i+1],lis_b[i],lis_y1r[i],self.K))
        lis_y1=lis_y1r
        lis_v=[]
        lis_E=[]
        g=9.807
        for i in range(len(lis_b)):
            A=lis_z[i]*lis_y1[i]**2+lis_b[i]*lis_y1[i]
            lis_v.append(self.Q/A)
            lis_E.append(lis_y1[i]+(lis_v[i]**2)/(2*g))    
        yf=max(lis_y1)
        L=int(L+1)
        L=4.7*b+1.65*self.zf*yf
        L=int(L+1)
        incx=L/self.nt
        x=[i*incx for i in range(self.nt+1)]
        if (show_results=="TRUE"):
            print("x[m]".center(12),"b[m]".center(12),"z".center(12),"HZ[m]".center(12),"y[m]".center(12),"v[m/s]".center(12),"E[m]".center(12))
            print("============================================================================================")
            for i in range(len(lis_v)):
                print(str(round(x[i],decimal_number)).center(12),str(round(lis_b[i],decimal_number)).center(12),
                      str(round(lis_z[i],decimal_number)).center(12),
                      str(round(lis_inch[i],decimal_number)).center(12),str(round(lis_y1[i],decimal_number)).center(12),
                      str(round(lis_v[i],decimal_number)).center(12),str(round(lis_E[i],decimal_number)).center(12))
        return [x,lis_b,lis_z,lis_inch,lis_y1,lis_v,lis_E,L,self.HZ]
    def show_transicion(self,lis_col=["cyan","black"]):
        lis_res=self.results(show_results="FALSE")
        lis_x=lis_res[0]
        lis_b=rev(lis_res[1])
        lis_z=rev(lis_res[2])
        lis_inch=rev(lis_res[3])
        lis_y1=rev(lis_res[4])
        lis_v=rev(lis_res[5])
        lis_E=rev(lis_res[6])
        L=lis_res[7]
        HZ=lis_res[8]
####################################
        bmin=min(lis_b)
        zmin=min(lis_z)
        ymin=min(lis_y1)
        bmax=max(lis_b)
        zmax=max(lis_z)
        ymax=max(lis_y1)
        desh=rev(lis_inch)
        points1 = np.array([
            [L,-0.5*bmin-zmin*ymin,ymin+self.HZ],
            [L,-0.5*bmin,0+self.HZ],
            [L,0.5*bmin,0+self.HZ],
            [L,0.5*bmin+zmin*ymin,ymin+self.HZ],
            [0,-0.5*bmax-zmax*ymax,ymax],
            [0,-0.5*bmax,0],
            [0,0.5*bmax, 0],
            [0,0.5*bmax+zmax*ymax,ymax]
            ])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        r = [-1,1]
        X, Y = np.meshgrid(r,r)
        # plot vertices
        ax.scatter3D(points1[:,0],points1[:,1],points1[:,2])
        verts1=[]
        for i in range(len(lis_b)):
            verts1.append([[lis_x[i],-0.5*lis_b[i]-lis_z[i]*lis_y1[i],lis_y1[i]+desh[i]],
                           [lis_x[i],-0.5*lis_b[i],0+desh[i]],
                           [lis_x[i],0.5*lis_b[i], 0+desh[i]],
                           [lis_x[i],0.5*lis_b[i]+lis_z[i]*lis_y1[i],lis_y1[i]+desh[i]]])
        ax.add_collection3d(Poly3DCollection(verts1,facecolors=lis_col[0],linewidths=1,edgecolors=lis_col[1], alpha=.25))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.title("Transición de entrada (contracción)")
        plt.show()

#Q=5#Caudal [m3/s]
#bf=6#Ancho de solera seccion f [m]
#bc=3#Ancho de solera seccion c [m]
#zf=2.5#Talud en la seccion f
#yc=2.050 #Tirante en la seccion c [m]
#HZ=0.1#Diferencia de cotas [m]
#K=0.1#Coeficiente de perdidas en la transicion
#nt=10
#rr=Transicion_AE(Q,bf,bc,zf,yc,HZ,K,nt)
#rr.results()
#rr.show_transicion()











