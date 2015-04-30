# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:00:29 2015

@author: Felipe
"""

import matplotlib.pyplot as plt
import math

def Fourier(d,A,k,Tparede_e,Tparede_i):
    Q=(k*A/d)*(Tparede_e-Tparede_i)
    return Q
    
def convec_entrada(h,A,Qentrada):
    dTi=Qentrada/(A*h)
    return dTi
    
def irradiação_leste(cos_leste):
    dT_leste=3.7*cos_leste
    return dT_leste

def irradiação_oeste(cos_oeste):
    dT_oeste=3.7*cos_oeste
    return dT_oeste

def irradiação_teto(cos_teto):
    dT_teto=6.2*cos_teto
    return dT_teto
    
def dU(n,R,Qentrada):
    dT=Qentrada/(2.5*n*R)
    return dT
    
def variação_em_Jaules(n,R,T):
    q=T*2.5*n*R
    return q
    
def energia_interna(n,R,K):
    q=2.5*n*R*K
    return q

def tempenergia(n,R,U):
    t=U/(n*R*2.5)
    return t
    
Posição_Sol=[0,0,0,0,0,0,5/180*math.pi,15/180*math.pi,30/180*math.pi,45/180*math.pi,60/180*math.pi,75/180*math.pi,90/180*math.pi,105/180*math.pi,120/180*math.pi,135/180*math.pi,150/180*math.pi,165/180*math.pi,175/180*math.pi,0,0,0,0,0]
Te=[0,-2,-5,-8,-10,-5,0,9,17,22,30,44,56,53,49,41,35,29,24,18,10,4,1,0]
T_Parede_oeste_Ext=[]*24
T_Parede_leste_Ext=[]*24
T_Teto=[]

Ti=[0]*24
Ti[0]=23

for i in range(0,24):
    if Posição_Sol[i]<math.pi/2:
        cos_teto1=math.sin(Posição_Sol[i])
        T_Teto.insert(i,cos_teto1)
        sen_parede_leste=math.cos(Posição_Sol[i])
        if sen_parede_leste!=1.0:           
            T_Parede_leste_Ext.insert(i,sen_parede_leste)
        if sen_parede_leste==1.0:
            T_Parede_leste_Ext.insert(i,0.0)
        T_Parede_oeste_Ext.insert(i,0.0)
    if Posição_Sol[i]>math.pi/2:
        cos_teto2=math.sin(Posição_Sol[i])
        T_Teto.insert(i,cos_teto2)
        sin_parede_oeste=math.cos(Posição_Sol[i])
        T_Parede_oeste_Ext.insert(i,sin_parede_oeste)
        T_Parede_leste_Ext.insert(i,0.0)
    if Posição_Sol[i]==math.pi/2:
        cos_teto3=math.sin(Posição_Sol[i])
        T_Teto.insert(i,cos_teto3)
        T_Parede_oeste_Ext.insert(i,0.0)
        T_Parede_leste_Ext.insert(i,0.0)

Tfinal_leste=[]        
Tfinal_oeste=[] 
Tfinal_norte=[]
Tfinal_sul=[]
Tfinal_teto=[]
for i in range(0,24):        
    dt_leste=irradiação_leste(T_Parede_leste_Ext[i])
    dt_oeste=-(irradiação_oeste(T_Parede_oeste_Ext[i]))
    dt_teto=irradiação_teto(T_Teto[i])
    a=Te[i]+dt_leste
    Tfinal_leste.append(a)
    b=Te[i]+dt_oeste
    Tfinal_oeste.append(b)
    c=Te[i]+dt_teto
    Tfinal_teto.append(c)
    Tfinal_norte.append(Te[i])
    Tfinal_sul.append(Te[i])

nAC=9.32 ##Numero de mols de ar que sai do ar condicionado
nAQ=58.7 #Numero de mols de ar que sai do aquecedor
KAC=280 #Temperatura em Kelvin do ar de saída do Ar Condicionado
KAQ=329 #Temperatura em Kelvin do ar de saída do Aquecedor
R=8.31 #constante universal dos gases
Mar=90000 #massa de ar em gramas
n=Mar/28.96 #3107
d=0.15 # espessura da parede
A=15 #área de cada parede
k=1.1 #k do fourier
Car=1.01 #calor especifico do ar
Ateto=25 #area do teto
h=[8,6,6,5,5,5,6,6,6,7,8,9,10,11,10,10,8,7,6,5,5,4,4,3] #h da equação de convecção
Atotal=85 #Área total

Quanto_deve_variar=[]
teste=[] #variação total de calor
teste1=[]#variação de temperatura
teste2=[]#variação de temperatura
teste3=[]#variação de temperatura
teste4=[]#variação de temperatura
teste5=[]#variação de temperatura
teste6=[]#variação de temperatura
Tempideal=[23]*24
jaules_segundo=[]
energiainternacasa=[]
energiainterna23=[]
deltaenergia=[]
deltaenergiasegAC=[]
deltaenergiasegAQ=[]
razaoAC23=[]
razaoAQ23=[]



for i in range(0,23):   
    Qteto=Fourier(d,Ateto,k,Tfinal_teto[i],Ti[i])
    teste2.append(Qteto)
    Qleste=Fourier(d,A,k,Tfinal_leste[i],Ti[i])
    teste3.append(Qleste)
    Qoeste=Fourier(d,A,k,Tfinal_oeste[i],Ti[i])
    teste4.append(Qoeste)
    Qsul=Fourier(d,A,k,Tfinal_sul[i],Ti[i])
    teste5.append(Qsul)
    Qnorte=Fourier(d,A,k,Tfinal_norte[i],Ti[i])
    teste6.append(Qnorte)
    Q_Total=Qteto+Qleste+Qoeste+Qsul+Qnorte
    teste.append(Q_Total)
    s=convec_entrada(h[i],Atotal,teste[i])
    teste1.append(s)
    jaules_segundo.append(teste[i]/3600)
    Ti[i+1]=Ti[i]+s
    aert=energia_interna(n,R,Ti[i]+273)
    energiainternacasa.append(aert)
    aer=energia_interna(n,R,23+273)
    energiainterna23.append(aer)
    deltaenergia.append(energiainternacasa[i]-energiainterna23[i])
    
for i in range(0,24):
    deltaenergiasegAC.append(energia_interna(nAC,R,Ti[i]+273))
    deltaenergiasegAQ.append(energia_interna(nAQ,R,Ti[i]+273))

    
    
jaules=[]
djaulesAC=[]
djaulesAQ=[]

for i in range(0,24):
    Quanto_deve_variar.append(Tempideal[i]-Ti[i])
    kk=variação_em_Jaules(n,R,Quanto_deve_variar[i])
    jaules.append(kk)
    
for i in range(0,24):
    if jaules[i]>0:
        pp=energia_interna(nAC,R,KAC)
        djaulesAC.append(pp)
    if jaules[i]<0:
        op=energia_interna(nAQ,R,KAQ)
        djaulesAQ.append(op)
contador0=0
contador1=0
contador2=0
contador3=0
contador4=0
contador5=0
contador6=0
contador7=0
contador8=0
contador9=0
contador10=0
contador11=0
contador12=0
contador13=0
contador14=0
contador15=0
contador16=0
contador17=0
contador18=0
contador19=0
contador20=0
contador21=0
contador22=0
contador23=0


tempfinaldps0=[23.00]
energiainternacasadps0=[0]*3600
energiainternacasadps0[0]=energiainternacasa[1]
deltaenergiadps0=[0]*3600
deltaenergiadps0[0]=deltaenergia[1]

for i in range(1,3600):
    if deltaenergiadps0[i-1]<0:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador0+=1
    if deltaenergiadps0[i-1]>=0:
        aaaa=deltaenergiasegAQ[0]
        bbbb=deltaenergiasegAQ[0]
    energiainternacasadps0[i]=energiainternacasadps0[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps0[i-1]+273))
    deltaenergiadps0[i]=deltaenergiadps0[i-1]+bbbb-deltaenergiasegAQ[0]
    rop=tempenergia(n,R,energiainternacasadps0[i])-273
    tempfinaldps0.append(rop)
energiainternacasadps1=[0]*3600
energiainternacasadps1[0]=energiainternacasadps0[3599]
deltaenergiadps1=[0]*3600
deltaenergiadps1[0]=deltaenergiadps0[3599]
tempfinaldps1=[tempfinaldps0[3599]]
  
for i in range(1,3600):
    if energiainternacasadps1[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador1+=1
    if energiainternacasadps1[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[1]
        bbbb=deltaenergiasegAQ[1]
    energiainternacasadps1[i]=energiainternacasadps1[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps1[i-1]+273))
    deltaenergiadps1[i]=deltaenergiadps1[i-1]+bbbb-deltaenergiasegAQ[1]
    rop=tempenergia(n,R,energiainternacasadps1[i])-273
    tempfinaldps1.append(rop)
energiainternacasadps2=[0]*3600
energiainternacasadps2[0]=energiainternacasadps1[3599]
deltaenergiadps2=[0]*3600
deltaenergiadps2[0]=deltaenergiadps1[3599]
tempfinaldps2=[tempfinaldps1[3599]]
    

for i in range(1,3600):
    if energiainternacasadps2[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador2+=1
    if energiainternacasadps2[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[2]
        bbbb=deltaenergiasegAQ[2]
    energiainternacasadps2[i]=energiainternacasadps2[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps2[i-1]+273))
    deltaenergiadps2[i]=deltaenergiadps2[i-1]+bbbb-deltaenergiasegAQ[2]
    rop=tempenergia(n,R,energiainternacasadps2[i])-273
    tempfinaldps2.append(rop)
energiainternacasadps3=[0]*3600
energiainternacasadps3[0]=energiainternacasadps2[3599]
deltaenergiadps3=[0]*3600
deltaenergiadps3[0]=deltaenergiadps2[3599]
tempfinaldps3=[tempfinaldps2[3599]]

for i in range(1,3600):
    if energiainternacasadps3[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador3+=1
    if energiainternacasadps3[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[3]
        bbbb=deltaenergiasegAQ[3]
    energiainternacasadps3[i]=energiainternacasadps3[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps3[i-1]+273))
    deltaenergiadps3[i]=deltaenergiadps3[i-1]+bbbb-deltaenergiasegAQ[3]
    rop=tempenergia(n,R,energiainternacasadps3[i])-273
    tempfinaldps3.append(rop)

energiainternacasadps4=[0]*3600
energiainternacasadps4[0]=energiainternacasadps3[3599]
deltaenergiadps4=[0]*3600
deltaenergiadps4[0]=deltaenergiadps3[3599]
tempfinaldps4=[tempfinaldps3[3599]]

for i in range(1,3600):
    if energiainternacasadps4[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador4+=1
    if energiainternacasadps4[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[4]
        bbbb=deltaenergiasegAQ[4]
    energiainternacasadps4[i]=energiainternacasadps4[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps4[i-1]+273))
    deltaenergiadps4[i]=deltaenergiadps4[i-1]+bbbb-deltaenergiasegAQ[4]
    rop=tempenergia(n,R,energiainternacasadps4[i])-273
    tempfinaldps4.append(rop)
energiainternacasadps5=[0]*3600
energiainternacasadps5[0]=energiainternacasadps4[3599]
deltaenergiadps5=[0]*3600
deltaenergiadps5[0]=deltaenergiadps4[3599]
tempfinaldps5=[tempfinaldps4[3599]]

for i in range(1,3600):
    if energiainternacasadps5[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador5+=1
    if energiainternacasadps5[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[5]
        bbbb=deltaenergiasegAQ[5]
    energiainternacasadps5[i]=energiainternacasadps5[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps5[i-1]+273))
    deltaenergiadps5[i]=deltaenergiadps5[i-1]+bbbb-deltaenergiasegAQ[5]
    rop=tempenergia(n,R,energiainternacasadps5[i])-273
    tempfinaldps5.append(rop)
energiainternacasadps6=[0]*3600
energiainternacasadps6[0]=energiainternacasadps5[3599]
deltaenergiadps6=[0]*3600
deltaenergiadps6[0]=deltaenergiadps5[3599]
tempfinaldps6=[tempfinaldps5[3599]]

for i in range(1,3600):
    if energiainternacasadps6[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador6+=1
    if energiainternacasadps6[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[6]
        bbbb=deltaenergiasegAQ[6]
    energiainternacasadps6[i]=energiainternacasadps6[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps6[i-1]+273))
    deltaenergiadps6[i]=deltaenergiadps6[i-1]+bbbb-deltaenergiasegAQ[6]
    rop=tempenergia(n,R,energiainternacasadps6[i])-273
    tempfinaldps6.append(rop)
    
    
energiainternacasadps7=[0]*3600
energiainternacasadps7[0]=energiainternacasadps6[3599]
deltaenergiadps7=[0]*3600
deltaenergiadps7[0]=deltaenergiadps6[3599]
tempfinaldps7=[tempfinaldps6[3599]]

for i in range(1,3600):
    if energiainternacasadps7[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador7+=1
    if energiainternacasadps7[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[7]
        bbbb=deltaenergiasegAQ[7]
    energiainternacasadps7[i]=energiainternacasadps7[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps7[i-1]+273))
    deltaenergiadps7[i]=deltaenergiadps7[i-1]+bbbb-deltaenergiasegAQ[7]
    rop=tempenergia(n,R,energiainternacasadps7[i])-273
    tempfinaldps7.append(rop)
    
energiainternacasadps8=[0]*3600
energiainternacasadps8[0]=energiainternacasadps7[3599]
deltaenergiadps8=[0]*3600
deltaenergiadps8[0]=deltaenergiadps7[3599]
tempfinaldps8=[tempfinaldps7[3599]]

for i in range(1,3600):
    if energiainternacasadps8[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador8+=1
    if energiainternacasadps8[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[8]
        bbbb=deltaenergiasegAQ[8]
    energiainternacasadps8[i]=energiainternacasadps8[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps8[i-1]+273))
    deltaenergiadps8[i]=deltaenergiadps8[i-1]+bbbb-deltaenergiasegAQ[8]
    rop=tempenergia(n,R,energiainternacasadps8[i])-273
    tempfinaldps8.append(rop)
    
energiainternacasadps9=[0]*3600
energiainternacasadps9[0]=energiainternacasadps8[3599]
deltaenergiadps9=[0]*3600
deltaenergiadps9[0]=deltaenergiadps8[3599]
tempfinaldps9=[tempfinaldps8[3599]]

for i in range(1,3600):
    if energiainternacasadps9[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador9+=1
    if energiainternacasadps9[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[9]
        bbbb=deltaenergiasegAQ[9]
    energiainternacasadps9[i]=energiainternacasadps9[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps9[i-1]+273))
    deltaenergiadps9[i]=deltaenergiadps9[i-1]+bbbb-deltaenergiasegAQ[9]
    rop=tempenergia(n,R,energiainternacasadps9[i])-273
    tempfinaldps9.append(rop)
    
energiainternacasadps10=[0]*3600
energiainternacasadps10[0]=energiainternacasadps9[3599]
deltaenergiadps10=[0]*3600
deltaenergiadps10[0]=deltaenergiadps9[3599]
tempfinaldps10=[tempfinaldps9[3599]]
listakipol=[]
for i in range(1,3600):
    kipol=(energia_interna(nAC,R,(tempfinaldps10[i-1]+273)))
    if energiainternacasadps10[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador10+=1
    if energiainternacasadps10[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[10]
        bbbb=0
    listakipol.append(kipol)
    energiainternacasadps10[i]=energiainternacasadps10[i-1]+aaaa-(kipol)
    deltaenergiadps10[i]=deltaenergiadps10[i-1]-bbbb+deltaenergiasegAC[10]
    rop=tempenergia(n,R,energiainternacasadps10[i])-273
    tempfinaldps10.append(rop)
    
energiainternacasadps11=[0]*3600
energiainternacasadps11[0]=energiainternacasadps10[3599]
deltaenergiadps11=[0]*3600
deltaenergiadps11[0]=deltaenergiadps10[3599]
tempfinaldps11=[tempfinaldps10[3599]]

for i in range(1,3600):
    kipol=(energia_interna(nAC,R,tempfinaldps11[i-1]+273))
    if energiainternacasadps11[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador11+=1
    if energiainternacasadps11[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[11]
        bbbb=deltaenergiasegAC[11]
    energiainternacasadps11[i]=energiainternacasadps11[i-1]+aaaa-kipol
    deltaenergiadps11[i]=deltaenergiadps11[i-1]+bbbb-deltaenergiasegAC[11]
    rop=tempenergia(n,R,energiainternacasadps11[i])-273
    tempfinaldps11.append(rop)
    
energiainternacasadps12=[0]*3600
energiainternacasadps12[0]=energiainternacasadps11[3599]
deltaenergiadps12=[0]*3600
deltaenergiadps12[0]=deltaenergiadps11[3599]
tempfinaldps12=[tempfinaldps11[3599]]

for i in range(1,3600):
    kipol=(energia_interna(nAC,R,tempfinaldps12[i-1]+273))
    if energiainternacasadps12[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador12+=1
    if energiainternacasadps12[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[12]
        bbbb=deltaenergiasegAC[12]
    energiainternacasadps12[i]=energiainternacasadps12[i-1]+aaaa-kipol
    rop=tempenergia(n,R,energiainternacasadps12[i])-273
    tempfinaldps12.append(rop)
    
energiainternacasadps13=[0]*3600
energiainternacasadps13[0]=energiainternacasadps12[3599]
deltaenergiadps13=[0]*3600
deltaenergiadps13[0]=deltaenergiadps12[3599]
tempfinaldps13=[tempfinaldps12[3599]]
testekipol=[]
for i in range(1,3600):
    kipol=(energia_interna(nAC,R,tempfinaldps13[i-1]+273))
    if energiainternacasadps13[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        testekipol.append(aaaa)
        bbbb=djaulesAC[1]
        contador13+=1
    if energiainternacasadps13[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[13]
        bbbb=deltaenergiasegAC[13]
    energiainternacasadps13[i]=energiainternacasadps13[i-1]+aaaa-kipol
    deltaenergiadps13[i]=deltaenergiadps13[i-1]-bbbb+deltaenergiasegAC[13]
    rop=tempenergia(n,R,energiainternacasadps13[i])-273
    tempfinaldps13.append(rop)
    
energiainternacasadps14=[0]*3600
energiainternacasadps14[0]=energiainternacasadps13[3599]
deltaenergiadps14=[0]*3600
deltaenergiadps14[0]=deltaenergiadps13[3599]
tempfinaldps14=[tempfinaldps13[3599]]

for i in range(1,3600):
    if energiainternacasadps14[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador14+=1
    if energiainternacasadps14[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[14]
        bbbb=deltaenergiasegAC[14]
    energiainternacasadps14[i]=energiainternacasadps14[i-1]+aaaa-(energia_interna(nAC,R,tempfinaldps14[i-1]+273))
    deltaenergiadps14[i]=deltaenergiadps14[i-1]-bbbb+deltaenergiasegAC[14]
    rop=tempenergia(n,R,energiainternacasadps14[i])-273
    tempfinaldps14.append(rop)

energiainternacasadps15=[0]*3600
energiainternacasadps15[0]=energiainternacasadps14[3599]
deltaenergiadps15=[0]*3600
deltaenergiadps15[0]=deltaenergiadps14[3599]
tempfinaldps15=[tempfinaldps14[3599]]

for i in range(1,3600):
    if energiainternacasadps15[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador15+=1
    if energiainternacasadps15[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[15]
        bbbb=deltaenergiasegAC[15]
    energiainternacasadps15[i]=energiainternacasadps15[i-1]+aaaa-(energia_interna(nAC,R,tempfinaldps15[i-1]+273))
    deltaenergiadps15[i]=deltaenergiadps15[i-1]+bbbb-deltaenergiasegAC[15]
    rop=tempenergia(n,R,energiainternacasadps15[i])-273
    tempfinaldps15.append(rop)
    
energiainternacasadps16=[0]*3600
energiainternacasadps16[0]=energiainternacasadps15[3599]
deltaenergiadps16=[0]*3600
deltaenergiadps16[0]=deltaenergiadps15[3599]
tempfinaldps16=[tempfinaldps15[3599]]

for i in range(1,3600):
    if energiainternacasadps16[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador16+=1
    if energiainternacasadps16[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[16]
        bbbb=deltaenergiasegAC[16]
    energiainternacasadps16[i]=energiainternacasadps16[i-1]+aaaa-(energia_interna(nAC,R,tempfinaldps16[i-1]+273))
    deltaenergiadps16[i]=deltaenergiadps16[i-1]-bbbb+deltaenergiasegAC[16]
    rop=tempenergia(n,R,energiainternacasadps16[i])-273
    tempfinaldps16.append(rop)
    
energiainternacasadps17=[0]*3600
energiainternacasadps17[0]=energiainternacasadps16[3599]
deltaenergiadps17=[0]*3600
deltaenergiadps17[0]=deltaenergiadps16[3599]
tempfinaldps17=[tempfinaldps16[3599]]

for i in range(1,3600):
    if energiainternacasadps17[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador17+=1
    if energiainternacasadps17[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[17]
        bbbb=deltaenergiasegAC[17]
    energiainternacasadps17[i]=energiainternacasadps17[i-1]+aaaa-(energia_interna(nAC,R,tempfinaldps17[i-1]+273))
    deltaenergiadps17[i]=deltaenergiadps17[i-1]-bbbb+deltaenergiasegAC[17]
    rop=tempenergia(n,R,energiainternacasadps7[i])-273
    tempfinaldps17.append(rop)
    
energiainternacasadps18=[0]*3600
energiainternacasadps18[0]=energiainternacasadps17[3599]
deltaenergiadps18=[0]*3600
deltaenergiadps18[0]=deltaenergiadps17[3599]
tempfinaldps18=[tempfinaldps17[3599]]

for i in range(1,3600):
    if energiainternacasadps18[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador18+=1
    if energiainternacasadps18[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[18]
        bbbb=deltaenergiasegAC[18]
    energiainternacasadps18[i]=energiainternacasadps18[i-1]+aaaa-(energia_interna(nAC,R,tempfinaldps18[i-1]+273))
    deltaenergiadps18[i]=deltaenergiadps18[i-1]-bbbb+deltaenergiasegAC[18]
    rop=tempenergia(n,R,energiainternacasadps18[i])-273
    tempfinaldps18.append(rop)
    
energiainternacasadps19=[0]*3600
energiainternacasadps19[0]=energiainternacasadps18[3599]
deltaenergiadps19=[0]*3600
deltaenergiadps19[0]=deltaenergiadps18[3599]
tempfinaldps19=[tempfinaldps18[3599]]
 
for i in range(1,3600):
    if energiainternacasadps19[i-1]>energiainterna23[1]:
        aaaa=djaulesAC[1]
        bbbb=djaulesAC[1]
        contador19+=1
    if energiainternacasadps19[i-1]<=energiainterna23[1]:
        aaaa=deltaenergiasegAC[19]
        bbbb=deltaenergiasegAC[19]
    energiainternacasadps19[i]=energiainternacasadps19[i-1]+aaaa-(energia_interna(nAC,R,tempfinaldps19[i-1]+273))
    deltaenergiadps19[i]=deltaenergiadps19[i-1]+bbbb-deltaenergiasegAC[19]
    rop=tempenergia(n,R,energiainternacasadps19[i])-273
    tempfinaldps19.append(rop)
    
energiainternacasadps20=[0]*3600
energiainternacasadps20[0]=energiainternacasadps19[3599]
deltaenergiadps20=[0]*3600
deltaenergiadps20[0]=deltaenergiadps19[3599]
tempfinaldps20=[tempfinaldps19[3599]]
     
for i in range(1,3600):
    if energiainternacasadps20[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador20+=1
    if energiainternacasadps20[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[20]
        bbbb=deltaenergiasegAQ[20]
    energiainternacasadps20[i]=energiainternacasadps20[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps20[i-1]+273))
    deltaenergiadps20[i]=deltaenergiadps20[i-1]+bbbb-deltaenergiasegAQ[20]
    rop=tempenergia(n,R,energiainternacasadps20[i])-273
    tempfinaldps20.append(rop)
    
energiainternacasadps21=[0]*3600
energiainternacasadps21[0]=energiainternacasadps20[3599]
deltaenergiadps21=[0]*3600
deltaenergiadps21[0]=deltaenergiadps20[3599]
tempfinaldps21=[tempfinaldps20[3599]]
  
for i in range(1,3600):
    if energiainternacasadps21[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador21+=1
    if energiainternacasadps21[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[21]
        bbbb=deltaenergiasegAQ[21]
    energiainternacasadps21[i]=energiainternacasadps21[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps21[i-1]+273))
    deltaenergiadps21[i]=deltaenergiadps21[i-1]+bbbb-deltaenergiasegAQ[21]
    rop=tempenergia(n,R,energiainternacasadps21[i])-273
    tempfinaldps21.append(rop)

energiainternacasadps22=[0]*3600
energiainternacasadps22[0]=energiainternacasadps21[3599]
deltaenergiadps22=[0]*3600
deltaenergiadps22[0]=deltaenergiadps21[3599]
tempfinaldps22=[tempfinaldps21[3599]]

for i in range(1,3600):
    if energiainternacasadps22[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador22+=1
    if energiainternacasadps22[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[22]
        bbbb=deltaenergiasegAQ[22]
    energiainternacasadps22[i]=energiainternacasadps22[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps22[i-1]+273))
    deltaenergiadps22[i]=deltaenergiadps22[i-1]+bbbb-deltaenergiasegAQ[22]
    rop=tempenergia(n,R,energiainternacasadps2[i])-273
    tempfinaldps22.append(rop)

energiainternacasadps23=[0]*3600
energiainternacasadps23[0]=energiainternacasadps22[3599]
deltaenergiadps23=[0]*3600
deltaenergiadps23[0]=deltaenergiadps22[3599]
tempfinaldps23=[tempfinaldps22[3599]]

for i in range(1,3600):
    if energiainternacasadps23[i-1]<energiainterna23[1]:
        aaaa=djaulesAQ[1]
        bbbb=djaulesAQ[1]
        contador23+=1
    if energiainternacasadps23[i-1]>=energiainterna23[1]:
        aaaa=deltaenergiasegAQ[23]
        bbbb=deltaenergiasegAQ[23]
    energiainternacasadps23[i]=energiainternacasadps23[i-1]+aaaa-(energia_interna(nAQ,R,tempfinaldps23[i-1]+273))
    deltaenergiadps23[i]=deltaenergiadps23[i-1]+bbbb-deltaenergiasegAQ[23]
    rop=tempenergia(n,R,energiainternacasadps3[i])-273
    tempfinaldps23.append(rop)
        
contadores=[]
contadores2=[]
contadores.append(contador0*0.44/1000)
contadores.append(contador1*0.44/1000)
contadores.append(contador2*0.44/1000)
contadores.append(contador3*0.44/1000)
contadores.append(contador4*0.44/1000)
contadores.append(contador5*0.44/1000)
contadores.append(contador6*0.44/1000)
contadores.append(contador7*0.44/1000)
contadores.append(contador8*0.44/1000)
contadores.append(contador9*0.44/1000)
contadores.append(contador10*-1.0/1000)        
contadores.append(contador11*-1.0/1000)
contadores.append(contador12*-1.0/1000)
contadores.append(contador13*-1.0/1000)
contadores.append(contador14*-1.0/1000)
contadores.append(contador15*-1.0/1000)
contadores.append(contador16*-1.0/1000)
contadores.append(contador17*-1.0/1000)
contadores.append(contador18*-1.0/1000)
contadores.append(contador19*-1.0/1000)
contadores.append(contador20*0.44/1000)
contadores.append(contador21*0.44/1000)
contadores.append(contador22*0.44/1000)
contadores.append(contador23*0.44/1000)

contadores2.append(contador0*-1.0/1000)
contadores2.append(contador1*-1.0/1000)
contadores2.append(contador2*-1.0/1000)
contadores2.append(contador3*-1.0/1000)
contadores2.append(contador4*-1.0/1000)
contadores2.append(contador5*-1.0/1000)
contadores2.append(contador6*-1.0/1000)
contadores2.append(contador7*-1.0/1000)
contadores2.append(contador8*-1.0/1000)
contadores2.append(contador9*-1.0/1000)
contadores2.append(contador10*0.5/1000)        
contadores2.append(contador11*0.5/1000)
contadores2.append(contador12*0.5/1000)
contadores2.append(contador13*0.5/1000)
contadores2.append(contador14*0.5/1000)
contadores2.append(contador15*0.5/1000)
contadores2.append(contador16*0.5/1000)
contadores2.append(contador17*0.5/1000)
contadores2.append(contador18*0.5/1000)
contadores2.append(contador19*0.5/1000)
contadores2.append(contador20*-1.0/1000)
contadores2.append(contador21*-1.0/1000)
contadores2.append(contador22*-1.0/1000)
contadores2.append(contador23*-1.0/1000)

soma=0.44/1000*(contador0+contador1+contador2+contador3+contador4+contador5+contador6+contador7+contador8+contador9)
soma1=0.5/1000*(contador10+contador11+contador12+contador13+contador14+contador15+contador16+contador17+contador18+contador19)
soma2=0.44/1000*(contador20+contador21+contador22+contador23)
somafinal=soma+soma1+soma2
Tempo=[0]*24
Tempo=[0]*24
for i in range(1,len(Tempo)):
    Tempo[i]=Tempo[i-1]+1

Tempos=[0]*3600
for u in range(1,len(Tempos)):
    Tempos[u]=Tempos[u-1]+1
    
    
Temposi=[0]*3599
for u in range(1,len(Temposi)):
    Temposi[u]=Temposi[u-1]+1
    
plt.plot(Tempo,Ti,'b',label='Temperatura Interna')
plt.plot(Tempo,Te,'r',label='Temperatura Externa')
plt.plot(Tempo,Tempideal,'green',linestyle='--',label='Temperatura Interna Ideal')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([0,23,-20,60])
plt.ylabel('Temperatura[°C]')
plt.xlabel('Horário')
plt.title('Temperaturas Interna e Externa sem o Termostato')
plt.show()


plt.plot(Tempo,Quanto_deve_variar,'blue',label="Delta Temperatura")
plt.plot(Tempo,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],'r',linestyle='--',label="Ponto de viragem AC-aquecedor/aquecedor-AC")
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.annotate('Troca Aquecedor-AC (10h)', xy=(9.9, 0), xytext=(0.2, -25),
            arrowprops=dict(facecolor='black', shrink=0.05),)
plt.annotate('Troca AC-Aquecedor (19h)', xy=(19.1, -0.5), xytext=(7, 25),
            arrowprops=dict(facecolor='black', shrink=0.05),)
plt.axis([0,23,-36,36])
plt.ylabel('Temperatura[°C]')
plt.xlabel('Horário')
plt.title('O quanto o termostato deve variar na temperatura para mante-la constante em 23°C')
plt.show()


plt.plot(Tempos,tempfinaldps8,'blue',label="Temperatura interna")
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([0,70,16,25])
plt.ylabel('Temperatura[°C]')
plt.xlabel('Horário[s]')
plt.title('Temperatura interna sob efeito do aquecedor na hora 1')
plt.show()
print("Nessa hora o gasto energético foi de ",contador8*0.5," Watts")
print("ou ", contador8*0.5/1000,'kW')

plt.plot(Tempo,contadores,'ro--',label="Aquecedor")
plt.plot(Tempo,contadores2,'bo--',label="Ar Condicionado")
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([0,23,0,1.4])
plt.ylabel('Energia[kW]')
plt.xlabel('Horário[h]')
plt.title('Energia gasta pelo termostato')
plt.show()
#
print("O gasto total de energia desse dia foi de",somafinal,"kW")
print("Convertendo para a moeda local, o gasto foi de",somafinal*1.14,"EGP")
print("Convertendo para dólares o gasto foi de",somafinal*0.15,"Dólares")