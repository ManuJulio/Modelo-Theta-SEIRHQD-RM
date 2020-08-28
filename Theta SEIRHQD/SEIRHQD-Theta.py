#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


#Fechas importantes:
#t_0= 01 de abril de 2020, para tener una base de datos mesurable
#24 de abril anuncio plan retorno a la nueva normalidad
#relajamiento medidas de control
#15 de mayo anuncio cuarentena total
#Ricardo Baeza-Yaetes propone hay 13 veces más activos que los registrados
# m: efectividad de las medidas [0,1]. 1 es muy inefectiva
# w: tasa de mortalidad, 
# se aumenta cuando el indice de saturación es superior a 1, en Chile es sobre 1000 contagiados
# yhd días de transicion de un hospitalizado a D, yhr a recuperado
# yid días de transicion de un infectado detectado a H
# teta: % de infectados identificados [0,1]. En 1 no hay casos sin detectar
# bid tasa de infeccion de un infectado detectado [0,0.5]
# be=Cebid, biu=Cubid, bhd=bhr=Chbid. bid tasa de infección de un infectado a otra etapa, Ce, Cu, Ch [0,1]
# bhd=bhr=tasa de infeccion de personal medico


# In[3]:


#25-02-2020	Entran los primeros dos contagiados
#29-02-2020	Entra el tercer y quinto contagiado al país
#03-03-2020	Se detecta al primer contagiado
#04-03-2020	Se detecta al segundo y tercer contagiado 
#05-03-2020	Se detecta al 4to contagio
#06-03-2020	Se detecta al 5to contagiado
#07-03-2020	6to contagiado entró con el 3er y 5to contagiado. 7mo contagio se detecta
#08-03-2020	Contagiada sale del país y deja una mujer contagiada de 83 años. Se llega a 10 contagios
#15-03-2020	Se suspenden las clases presenciales
#18-03-2020	Se declara cierre de fronteras desde las 00:00 del 19 de marzo por 90 días
#	Hay 8 hospitalizados, 3 graves.
#	Se modifica requisito de caso sospechoso
#	La Reina anuncia cuarentena para el 20 de marzo
#	Se decreta cierre de centros comerciales para el 19 de marzo
#19-03-2020	Se confirman 103 nuevos casos y 342 acumulados
#	247 se encuentran en RM
#	Se declara aduana sanitaria desde el 20 de marzo para evitar contagio en algunas regiones
#	19 hospitalizados, 5 recuperados
#20-03-2020	Walmart adelanta el cierro de locales a las 20:00
#	56 alcaldes solicitan cuarentena nacional al presidente
#	Maipu declara cuarentena total en la comuna
#	La Reina, Las Condes y Vitacura anuncian cuarentena preventiba para el 21 de marzo
#	Se anuncia cierre de cines, teatros, pubs, bares, discotecas, restorantes y centros deportivos en el país
#21-03-2020	33 pacientes hospitalizados y 8 recuperados
#	Primer fallecido, mujer de 83 años de Renca
#23-03-2020	Segunda fallecida por Covid, de Maipú
#24-03-2020	Se fija el precio del examen en 25.000 pesos
#	Se llama a los centros médicos a postergar cirujías de carácter electivo
#26-06-2020	Cuarentena total preventiva para Vitacura, Las Condes, Lo Barnechea, Providencia, Santiago, Ñuñoa e Independencia
#	Se exime a empleadores de pagar a los trabajadores por la cuarentena
#29-03-2020	Primera disputa en el conteo de fallecidos
#01-04-2020	Inician controles sanitarios en los terminales de buses del país
#	293 nuevos casos, 3031 acumulados, 1521 en RM. 232 recuperados, 142 hospitalizados
#	Gobierno se querella con dos alcades de Ñuble por cerrar el acceso
#	
#02-04-2020	Finaliza cuarentena en Independencia 22:00
#03-04-2020	 Se declara cordon sanitario para RM y Concepción para el 9 a 12 de abril (semana santa)
#04-04-2020	424 nuevos casos, 4161 acumulados, 528 recuperados, 1957 en RM, 
#	280 en UCI, 225 con ventilador, 38 en estado critico. 27 fallecidos
#05-04-2020	Médico de 61 años se contagia
#06-04-2020	Alcaldesa de Antofagasta solicita al presidente declarar cuarentena en la región
#	Las Condes decreta ordenanza de uso obligatorio de mascarillas en la comuna
#	
#07-04-2020	301 nuevos contagios, 5116 acumulados.
#08-04-2020	Gobierno decreta uso obligatorio de mascarilla en el transporte público
#09-04-2020	Parte de Puente Alto entra en cuarentena (cuarentenas parciales geograficas 22:00)
#11-04-2020	Incremento sustancial de nuevos casos en Puente Alto
#13-04-2020	Finaliza cuarentena en Vitacura, Lo Barnechea, Providencia y Sec Norte de Ñuñoa y Stgo
#	312 nuevos contagios, 7525 acumulados, 308 recuperados y dos nuevos fallecidos


# In[ ]:





# In[32]:


def SEIHQRD(y,t,N,ye,yinf):
    S,E,I,Iu,Idu,Hr,Hd,Q,Rd,Ru,Du,D=y
    me=0.58889915 #medida de control aplicada a los expuestos en t0 (1 de abril de 2020)
    miu=0.418343253 #medida de control aplicada a los infectados no detectados en t0
    mid=0.41370334 #medida de control aplicada a los infectados detectados en t0
    mhr=0.2994201 #medida de control aplicada a los hospitalizados que se recuperaran en t0
    mhd=0.29990032 #medida de control aplicada a los hospitalizados que falleceran en t0
    bid=0.302230547 #tasa de contagio por estar en contacto con alguien infectado
    ce=0.903887292
    cu=0.871018708
    ch=0.139375658
    be=ce*bid #tasa de contagio por estar en contacto con alguien expuesto
    biu=cu*bid #tasa de contagio por estar en contacto con alguien asintomatico
    bhd=ch*bid #tasa de contagio por estar en contacto con un hospitalizado
    bhr=bhd
    yhd=0.037475687 #días de transición (días**-1) de un hospitalizado a estado D
    yhr=0.078839256 #días de transición (días**-1) de un hospitalizado a estado Q
    yi=1 #días de transición (días**-1) de un infectado a estado a Iu,Q,Hr,Hd
    yiu=1/9 #días de transición (días**-1) de un infectado no detectado a estado Ru
    w=0.019 #Tasa de fatalidad del virus
    wu=w #Tasa de fatalidad para no detectados
    yidu=1/10 #días de transición (días**-1) de un infectado no detectado al estado D
    teta=0.35837 #Proporción de infectados detectados
    yq=1/14 #días de transición (días**-1) de alguien en cuarentena al estado R (días impuestos por autoridad)
    p=0.028 #probabilidad de un infectado de ser hospitalizado
    if t>28.01: #a partir del 28 de abril se aplica relajo de medidas por llamado a la nueva normalidad
        teta=0.7229 #Proporción de infectados detectados cambia en función de saturación sistema atención
        me=0.652410404
        miu=0.432222681
        mid=0.403307236
        mhr=0.299865256
        mhd=0.300033437
        ce=0.873104031
        cu=0.869448088
        ch=0.137222316
        bid=0.272203775
        yhd=1/14 #días de un hospitalizado a D disminuyen a medida que se satura el sistema
        be=ce*bid
        biu=cu*bid
        bhd=ch*bid
        bhr=bhd
        w=0.022 #tasa de fatalidad aumenta con saturación del sistema
    if t>41.01: #11 de mayo se reaplican medidas tras aumento de contagios en abril
        teta=0.7080
        me=0.693554262
        miu=0.482486993
        mid=0.466126264
        mhr=0.199387429
        mhd=0.199931628
        ce=0.898357762
        cu=0.870867463
        ch=0.139281752
        bid=0.238821633
        be=ce*bid
        biu=cu*bid
        bhd=ch*bid
        bhr=bhd
        yhd=1/12
        w=0.02622
    if t>55.01:
        yhd=1/7
    if t>71.01:
        teta=0.8139
        me=0.587470322
        miu=0.332486993
        mid=0.316126264
        mhr=0.300223233
        mhd=0.300042725
        ce=0.898357762
        cu=0.870867463
        ch=0.139281752
        bid=0.242728381
        be=ce*bid
        biu=cu*bid
        bhd=ch*bid
        bhr=bhd
        w=0.02457
    if t>95: #14 de julio caida de fatalidad en espera de nuevos datos
        teta=0.8139
        me=0.507470322
        miu=0.252486993
        mid=0.236126264
        mhr=0.230223233
        mhd=0.230042725
        ce=0.898357762
        cu=0.870867463
        ch=0.139281752
        bid=0.242728381
        be=ce*bid
        biu=cu*bid
        bhd=ch*bid
        bhr=bhd
        w=0.009
    dS_dt=-(S/N)*(me*be*E+miu*biu*Iu+mid*bid*I+mhr*bhr*Hr+mhd*bhd*Hd)
    dE_dt=(S/N)*(me*be*E+miu*biu*Iu+mid*bid*I+mhr*bhr*Hr+mhd*bhd*Hd)-ye*E
    dI_dt=ye*E-yi*I
    dIu_dt=(1-teta-wu)*yi*I-yiu*Iu
    dIdu_dt=wu*yi*I-yidu*Idu
    dHr_dt=p*(teta-w)*yi*I-yhr*Hr
    dHd_dt=w*yi*I-yhd*Hd
    dQ_dt=(1-p)*(teta-w)*yi*I+yhr*Hr-yq*Q
    dRd_dt=yq*Q
    dRu_dt=yiu*Iu
    dDu_dt=yidu*Idu
    dD_dt=yhd*Hd
    return(dS_dt,dE_dt,dI_dt,dIu_dt,dIdu_dt,dHr_dt,dHd_dt,dQ_dt,dRd_dt,dRu_dt,dDu_dt,dD_dt)


# In[33]:


S=8119853
E=1828
I=320
Iu=594
Hr=81
Hd=8
Rd=234
Ru=200
D=5
N=8125072
Q=1930
Idu=12
Du=7
ye=0.1818
yinf=0.0714
t=np.linspace(0,275,6500)
ans=scipy.integrate.odeint(SEIHQRD,[S,E,I,Iu,Idu,Hr,Hd,Q,Rd,Ru,Du,D],t,args=(N,ye,yinf))
ans=np.array(ans)


# In[34]:


plt.figure(figsize=[15,6])
plt.plot(t,ans[:,6],label="Hd(t)")
plt.plot(t,ans[:,11],label="D(t)")
plt.plot(118, 7371.72,'ro',label='Fallecidos hoy')
plt.plot(118, 120.94,'go',label='Hospitalizados que falleceran')
plt.axvline(x=119, color='m',label='Paso 1')
plt.legend()
plt.xlabel("Tiempo")
plt.show()


# In[35]:


plt.figure(figsize=[15,6])
plt.plot(t,ans[:,2],label="I(t)")
plt.plot(t,ans[:,3],label="Iu(t)")
plt.axvline(x=119, color='m',label='Paso 2')
plt.legend()
plt.xlabel("Tiempo")
plt.show()


# In[26]:


plt.figure(figsize=[15,6])

plt.plot(t,ans[:,2],label="I(t)")

plt.axvline(x=119, color='m',label='Paso 2')
plt.legend()
plt.xlabel("Tiempo")
plt.show()


# In[22]:


plt.figure(figsize=[15,6])
plt.plot(t,ans[:,2],label="I(t)")
plt.axvline(x=119, color='m',label='Paso 2')
plt.legend()
plt.xlabel("Tiempo")
plt.show()


# In[21]:


plt.figure(figsize=[15,6])
plt.plot(t,ans[:,5],label="Hr(t)")
plt.plot(t,ans[:,6],label="Hd(t)")
plt.legend()
plt.xlabel("Tiempo")
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




