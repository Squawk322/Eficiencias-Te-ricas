import sys,time
import pandas as pd
import math as mt
import numpy as np
from PyQt5 import uic,QtWidgets,QtGui
from os import listdir, getcwd,path
import pyqtgraph as pg
from scipy import interpolate 
#%%
form_class = uic.loadUiType("puntual.ui")[0]
form_class_2= uic.loadUiType("seleccion.ui")[0]

def legencoef(n):
    p0 = np.array([1])
    p1 = np.array([1,0])
    if n==0:
        p=p0
    elif n==1:
        p=p1
    else:
        for i in range(2,n+1):
            pn = ((2*i-1)*np.append(p1,0)-(i-1)*np.append([0,0],p0))/i
            p0=p1
            p1=pn
        p=pn
    xl=np.sort(np.roots(p))
    w=np.zeros(n)
    for i in range(n):
        yl=np.zeros(n)
        yl[i]=1
        p=np.polyfit(xl,yl,n-1)
        P=np.polyint(p)
        w[i]=np.polyval(P,1)-np.polyval(P,-1)
    return xl,w

def read_coef(imput_file,):
    return pd.read_csv(imput_file, sep='\s+',header=None, names=['Energias','Coe_tot','Coe_fot'] )

def main(co,d,r,l,xv,denv,den,met,orden):
    coeal=read_coef('coeal.txt')
    coenai=read_coef('coenai2.txt')#para tener hast 20 mv
    Eal=np.array(coeal['Energias'])
    Enai=np.array(coenai['Energias'])
    ual=np.array(coeal['Coe_tot'])
    unai=np.array(coenai['Coe_fot'])
    E=np.unique(np.concatenate((Eal,Enai)))
    e = np.zeros([len(d),len(E)])
    eang = 1.0/(2*np.pi)
    uvent=interpolate.interp1d(Eal,ual,fill_value='extrapolate')(E)  #uvent=np.interp(E,Eal,ual)
    udetec=interpolate.interp1d(Enai,unai,fill_value='extrapolate')(E)  #udetec=np.interp(E,Enai,unai)
#    1: Trapecio   2:Simpson   3:Gauss-Legendre   4:MC
    if met==1:
        n=128
        #print('trapecio')
        for dn,dis in enumerate(d):
            for En,Et in enumerate(E):
                uventana=uvent[En]
                udetector=udetec[En]
                ec=np.zeros(n+1)
                ar=np.arange(n+1)
                dist=dis+xv
                if co>=0 and co<=r:
                    hfi=mt.pi/n
                    fi=ar*hfi
                    for finum,fival in enumerate(fi):
                        g=(mt.sqrt(r**2-(co*mt.sin(fival))**2)-co*mt.cos(fival))
                        a=mt.atan(g/(dist+l))
                        b=mt.atan(g/dist)
                        e1=0
                        e2=0
                        if 0<a:
                            h1=a/n
                            te=ar*h1
                            x=l/np.cos(te)
                            f1=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e1=h1*(f1.sum()-(f1[0]+f1[-1])/2)
                        if a<b:
                            h2=(b-a)/n
                            te=a+ar*h2
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=h2*(f2.sum()-(f2[0]+f2[-1])/2)
                        ec[finum]=e1+e2
                else:
                    hfi=mt.asin(r/co)/n
                    fi=ar*hfi
                    for finum,fival in enumerate(fi):
                        g=(co*mt.cos(fival)+mt.sqrt(r**2-(co*mt.sin(fival))**2))
                        g2=(co*mt.cos(fival)-mt.sqrt(r**2-(co*mt.sin(fival))**2))
                        b=mt.atan(g/(dist+l))
                        if dist==0:
                            a=mt.pi/2
                            c=a
                        else:
                            a=mt.atan(g2/dist)
                            c=mt.atan(g/dist)
                        e1=0
                        e2=0
                        if a<b:
                            h1=(b-a)/n
                            te=a+ar*h1
                            x=l/np.cos(te)
                            f1=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e1=h1*(f1.sum()-(f1[0]+f1[-1])/2)
                        if b<c and a<b:
                            h2=(c-b)/n
                            te=b+ar*h2
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=h2*(f2.sum()-(f2[0]+f2[-1])/2)
                        if a>b and a<c:
                            h2=(c-a)/n
                            te=a+ar*h2
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=h2*(f2.sum()-(f2[0]+f2[-1])/2)
                        ec[finum]=e1+e2
                e12=hfi*(ec.sum()-(ec[0]+ec[-1])/2)
                e[dn,En]=e12*eang*100    
    elif met==2:
        n=32
        #print('simpson')
        for dn,dis in enumerate(d):
            for En,Et in enumerate(E):
                uventana=uvent[En]
                udetector=udetec[En]
                ec=np.zeros(n+1)
                ar=np.arange(n+1)
                dist=dis+xv
                if co>=0 and co<=r:
                    hfi=mt.pi/n
                    fi=ar*hfi
                    for finum,fival in enumerate(fi):
                        g=(mt.sqrt(r**2-(co*mt.sin(fival))**2)-co*mt.cos(fival))
                        a=mt.atan(g/(dist+l))
                        b=mt.atan(g/dist)
                        e1=0
                        e2=0
                        if 0<a:
                            h1=a/n
                            te=ar*h1
                            x=l/np.cos(te)
                            f1=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e1=(2*f1[::2].sum()+4*f1[1::2].sum()-f1[0]-f1[-1])*h1/3
                        if a<b:
                            h2=(b-a)/n
                            te=a+ar*h2
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=(2*f2[::2].sum()+4*f2[1::2].sum()-f2[0]-f2[-1])*h2/3
                        ec[finum]=e1+e2
                else:
                    hfi=mt.asin(r/co)/n
                    fi=ar*hfi
                    for finum,fival in enumerate(fi):
                        g=(co*mt.cos(fival)+mt.sqrt(r**2-(co*mt.sin(fival))**2))
                        g2=(co*mt.cos(fival)-mt.sqrt(r**2-(co*mt.sin(fival))**2))
                        b=mt.atan(g/(dist+l))
                        if dist==0:
                            a=mt.pi/2
                            c=a
                        else:
                            a=mt.atan(g2/dist)
                            c=mt.atan(g/dist)
                        e1=0
                        e2=0
                        if a<b:
                            h1=(b-a)/n
                            te=a+ar*h1
                            x=l/np.cos(te)
                            f1=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e1=(2*f1[::2].sum()+4*f1[1::2].sum()-f1[0]-f1[-1])*h1/3
                        if b<c and a<b:
                            h2=(c-b)/n
                            te=b+ar*h2
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=(2*f2[::2].sum()+4*f2[1::2].sum()-f2[0]-f2[-1])*h2/3
                        if a>b and a<c:
                            h2=(c-a)/n
                            te=a+ar*h2
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=(2*f2[::2].sum()+4*f2[1::2].sum()-f2[0]-f2[-1])*h2/3
                        ec[finum]=e1+e2
                e12=(2*ec[::2].sum()+4*ec[1::2].sum()-ec[0]-ec[-1])*hfi/3
                e[dn,En]=e12*eang*100             
    elif met==3:
        xl,w=legencoef(orden)
        #print('gausslegendre')
        for dn,dis in enumerate(d):
            for En,Et in enumerate(E):
                uventana=uvent[En]
                udetector=udetec[En]
                ec=np.zeros(orden)
                dist=dis+xv
                if co>=0 and co<=r:
                    fi=0.5*(mt.pi*xl+mt.pi)
                    for finum,fival in enumerate(fi):
                        g=(mt.sqrt(r**2-(co*mt.sin(fival))**2)-co*mt.cos(fival))
                        a=mt.atan(g/(dist+l))
                        b=mt.atan(g/dist)
                        e1=0
                        e2=0
                        if 0<a:
                            te=0.5*(a*xl+a)
                            x=l/np.cos(te)
                            f1=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e1=sum(w*f1)*a/2
                        if a<b:
                            te=0.5*((b-a)*xl+a+b)
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=sum(w*f2)*(b-a)/2
                        ec[finum]=e1+e2
                else:
                    fi=0.5*(mt.asin(r/co)*xl+mt.asin(r/co))
                    for finum,fival in enumerate(fi):
                        g=(co*mt.cos(fival)+mt.sqrt(r**2-(co*mt.sin(fival))**2))
                        g2=(co*mt.cos(fival)-mt.sqrt(r**2-(co*mt.sin(fival))**2))
                        b=mt.atan(g/(dist+l))
                        if dist==0:
                            a=mt.pi/2
                            c=a
                        else:
                            a=mt.atan(g2/dist)
                            c=mt.atan(g/dist)
                        e1=0
                        e2=0
                        if a<b:
                            te=0.5*((b-a)*xl+a+b)
                            x=l/np.cos(te)
                            f1=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e1=sum(w*f1)*(b-a)/2
                        if b<c and a<b:
                            te=0.5*((c-b)*xl+b+c)
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=sum(w*f2)*(c-b)/2
                        if a>b and a<c:
                            te=0.5*((c-a)*xl+a+c)  
                            x=g/np.sin(te)-dist/np.cos(te)
                            f2=(1-np.exp(-udetector*den*x))*np.sin(te)*np.exp(-uventana*denv*xv/np.cos(te))
                            e2=sum(w*f2)*(c-a)/2
                        ec[finum]=e1+e2
                e12=sum(w*ec)*mt.pi/2
                e[dn,En]=e12*eang*100        
    elif met==4:
        #print('montecarlo')
        for dn,dist in enumerate(d):
            for En,Et in enumerate(E):
                ndet=nsup=ning=0 #detectados, superficie ventana (geometrica), ingresan al cristal
                uventana=uvent[En]
                udetector=udetec[En]
                div=0
                k=1
                num=100000
                if orden>num:
                    div=int(orden/num)
                    while k<=div:
                        one=np.ones(num)
                        x=co*one
                        y=np.zeros(num)
                        te=np.arccos(np.random.rand(num))
                        fi=2*mt.pi*np.random.rand(num)
                        x=x+dist*np.tan(te)*np.cos(fi);y=y+dist*np.tan(te)*np.sin(fi)
                        d1=r*r
                        remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                        nsup=nsup+len(x)
                        ve=-np.log(1-np.random.rand(len(x)))/(uventana*denv)
                        remove = [a for a, val in enumerate(ve*np.cos(te)) if val<=xv]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                        x=x+xv*np.tan(te)*np.cos(fi);y=y+xv*np.tan(te)*np.sin(fi)
                        remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                        ning=ning+len(x)
                        na=-np.log(1-np.random.rand(len(x)))/(udetector*den)
                        remove = [a for a, val in enumerate(na*np.cos(te)) if val>l]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove);na=np.delete(na,remove)
                        x=x+na*np.tan(te)*np.cos(fi);y=y+na*np.tan(te)*np.sin(fi)
                        remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                        ndet=ndet+len(x)
                    if orden-div*num>0:
                        num=orden-div*num
                        one=np.ones(num)
                        x=co*one
                        y=np.zeros(num)
                        te=np.arccos(np.random.rand(num))
                        fi=2*mt.pi*np.random.rand(num)
                        x=x+dist*np.tan(te)*np.cos(fi);y=y+dist*np.tan(te)*np.sin(fi)
                        d1=r*r
                        remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                        nsup=nsup+len(x)
                        ve=-np.log(1-np.random.rand(len(x)))/(uventana*denv)
                        remove = [a for a, val in enumerate(ve*np.cos(te)) if val<=xv]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                        x=x+xv*np.tan(te)*np.cos(fi);y=y+xv*np.tan(te)*np.sin(fi)
                        remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                        ning=ning+len(x)
                        na=-np.log(1-np.random.rand(len(x)))/(udetector*den)
                        remove = [a for a, val in enumerate(na*np.cos(te)) if val>l]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove);na=np.delete(na,remove)
                        x=x+na*np.tan(te)*np.cos(fi);y=y+na*np.tan(te)*np.sin(fi)
                        remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                        x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                        ndet=ndet+len(x)
                else:
                    num=orden
                    one=np.ones(num)
                    x=co*one
                    y=np.zeros(num)
                    te=np.arccos(np.random.rand(num))
                    fi=2*mt.pi*np.random.rand(num)
                    x=x+dist*np.tan(te)*np.cos(fi);y=y+dist*np.tan(te)*np.sin(fi)
                    d1=r*r
                    remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                    x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                    nsup=nsup+len(x)
                    ve=-np.log(1-np.random.rand(len(x)))/(uventana*denv)
                    remove = [a for a, val in enumerate(ve*np.cos(te)) if val<=xv]
                    x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                    x=x+xv*np.tan(te)*np.cos(fi);y=y+xv*np.tan(te)*np.sin(fi)
                    remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                    x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                    ning=ning+len(x)
                    na=-np.log(1-np.random.rand(len(x)))/(udetector*den)
                    remove = [a for a, val in enumerate(na*np.cos(te)) if val>l]
                    x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove);na=np.delete(na,remove)
                    x=x+na*np.tan(te)*np.cos(fi);y=y+na*np.tan(te)*np.sin(fi)
                    remove = [a for a, val in enumerate(x*x+y*y) if val>=d1]
                    x=np.delete(x,remove);y=np.delete(y,remove);te=np.delete(te,remove);fi=np.delete(fi,remove)
                    ndet=ndet+len(x)
                e[dn,En]=ndet*50/orden
    return e,E

class Seleccion(QtWidgets.QDialog,form_class_2):
     def __init__(self, parent=None):
          QtWidgets.QDialog.__init__(self, parent)
          self.setupUi(self)
          self.calc.clicked.connect(self.calc_clicked)
          self.save.clicked.connect(self.save_clicked)
          self.cancel.clicked.connect(self.cancel_clicked)
          self.onlynum=QtGui.QDoubleValidator()
          self.corrimiento.setValidator(self.onlynum)
          self.dini.setValidator(self.onlynum)
          self.dfin.setValidator(self.onlynum)
          self.radio.setValidator(self.onlynum)
          self.largo.setValidator(self.onlynum)
          self.grosorventana.setValidator(self.onlynum)
          self.denventana.setValidator(self.onlynum)
          self.dendetector.setValidator(self.onlynum)
          self.trap.clicked.connect(self.trap_clicked)
          self.simp.clicked.connect(self.simp_clicked)
          self.gauss.clicked.connect(self.gauss_clicked)
          self.mc.clicked.connect(self.mc_clicked)
          self.trap_clicked()
        
     def trap_clicked(self):
         self.metodo=1
         self.orden.setEnabled(False)
         self.trap.setChecked(True)
         self.intentos.setEnabled(False)
         
     def simp_clicked(self):
         self.metodo=2
         self.orden.setEnabled(False)
         self.intentos.setEnabled(False)
 
     def gauss_clicked(self):
         self.metodo=3
         self.orden.setEnabled(True)
         self.intentos.setEnabled(False)

     def mc_clicked(self):
         self.metodo=4
         self.orden.setEnabled(False)
         self.intentos.setEnabled(True)
                  
     def calc_clicked(self):
         co=abs(float(self.corrimiento.text()))
         inter=self.dint.value()
         d=np.arange(float(self.dini.text()),float(self.dfin.text())+inter,inter)
         self.d2=np.array([np.append([0],d)])
         r=float(self.radio.text())
         l=float(self.largo.text())
         xv=float(self.grosorventana.text())
         denv=float(self.denventana.text())
         den=float(self.dendetector.text())
         orden=1
         if self.gauss.isChecked():
             orden=self.orden.value()
         if self.mc.isChecked():
             orden=self.intentos.value()
         tiempo=time.time()
         self.a,self.b=main(co,d,r,l,xv,denv,den,self.metodo,orden)
         self.tiempo.setText("{} segundos".format(round(time.time() - tiempo,2)))
         self.noti.setText('Cálculos Completados')
         self.save.setEnabled(True)
    
     def save_clicked(self):
        file=self.filename.text()+'.txt'
        if file=='.txt':
            self.noname.setText('Ingrese un nombre de archivo')
        else:
            A=np.concatenate((self.d2,np.concatenate((np.array([self.b]),self.a)).T))
            np.savetxt('database/'+file, A)
            self.cancel_clicked()
     
     def cancel_clicked(self):
           self.noti.setText('No hay datos para guardar')
           self.tiempo.clear()
           self.filename.clear()
           self.noname.clear()
           self.save.setEnabled(False)
           self.trap_clicked()
           self.close()
                       
class MyWindowClass(QtWidgets.QMainWindow,form_class):
     def __init__(self, parent=None):
          QtWidgets.QMainWindow.__init__(self, parent)
          self.setupUi(self)
          self.seleccion=Seleccion()
          self.graph.clicked.connect(self.graph_clicked)
          self.addfile.clicked.connect(self.addfile_clicked)
          self.openfile.clicked.connect(self.openfile_clicked)
          self.refresh.clicked.connect(self.refresh_clicked)
          self.refresh_clicked()
          self.grafico.setTitle('EFICIENCIA vs ENERGIA')
          self.grafico.setLabel('bottom', 'Energia', units='eV')
          self.grafico.setLabel('left', 'Eficiencia (%)')
          self.grafico.setLogMode(x=False,y=False)
          self.xnor.clicked.connect(self.xnor_clicked)
          self.xlog.clicked.connect(self.xlog_clicked)
          self.ynor.clicked.connect(self.ynor_clicked)
          self.ylog.clicked.connect(self.ylog_clicked)
          self.onlynum=QtGui.QDoubleValidator()
          self.Eini.setValidator(self.onlynum)
          self.Efin.setValidator(self.onlynum)
          #self.showtab.clicked.connect(self.showtab_clicked)
     
     def xnor_clicked(self):
         self.grafico.setLogMode(x=False)
    
     def xlog_clicked(self):
         self.grafico.setLogMode(x=True)
    
     def ynor_clicked(self):
         self.grafico.setLogMode(y=False)
    
     def ylog_clicked(self):
         self.grafico.setLogMode(y=True)             
           
     def addfile_clicked(self):
         self.seleccion.exec_()
         self.refresh_clicked() 
         
     def refresh_clicked(self):
         self.ruta=getcwd()+'\\database'
         nombres=[path.splitext(x)[0] for x in listdir(self.ruta) if x.endswith('.txt')]
         self.lista.clear()
         if len(nombres)!=0:
             self.lista.setEnabled(True) 
             self.openfile.setEnabled(True)
             self.lista.addItem('Seleccione uno de los Archivos de Datos')
             for file2 in nombres:
                 self.lista.addItem(file2)
         else:
             self.lista.addItem('No hay Archivos de Datos Almacenados')
             self.lista.setEnabled(False)
             self.openfile.setEnabled(False)
    
     def openfile_clicked(self):
         self.grafico.clear()
         self.graph.setEnabled(False)
         filename2=self.lista.currentText()+'.txt'
         if filename2!='Seleccione uno de los Archivos de Datos.txt' and filename2!='No hay Archivos de Datos Almacenados.txt':
             data=np.loadtxt(self.ruta+'\\'+filename2)
             #dist=data[0]
             data=np.delete(data,0,0).T
             self.Ener=data[0]
             self.efis=np.delete(data,0,0)
             self.noti.setText('Archivo '+filename2+' Cargado Correctamente')
             self.graph.setEnabled(True)
          
     def graph_clicked(self):
         self.grafico.clear()
         inte=self.Eint.value()
         Eped=np.arange(float(self.Eini.text()),float(self.Efin.text())+inte,inte)
         #Eped=np.unique(np.concatenate((Eped,self.Ener)))
         lend=len(self.efis)
         lenE=len(Eped)
         c=np.zeros((lend,lenE))
         for dn,ef in enumerate(self.efis):
             c[dn]=np.interp(Eped,self.Ener,ef)
             self.grafico.plot(Eped*1000000, c[dn], pen=(dn,3))
         self.grafico.showGrid(x=True, y=True)

     #def showtab_clicked(self):
                 
         
         

#%% 
if __name__=="__main__":
    app = QtWidgets.QApplication(sys.argv)
    MyWindow = MyWindowClass(None)
    MyWindow.show()
    app.exec_()
