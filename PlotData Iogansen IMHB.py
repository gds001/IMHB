import numpy as np
import matplotlib.pyplot as plt
import json
import sklearn.linear_model as lm
import scipy.optimize as opt
def RMSD(e,f):
    return np.sqrt(np.sum((e-f)**2)/len(e))
def MSE(e,f):
    return np.average(e-f)
def RMSDf(params,hb,m,e,func):
    return np.sqrt(np.sum((e-func(hb,m,*params))**2)/len(e))
def Iogansen(hb,m,a,b,c):
    return a*np.sqrt(hb+c*m+b)-a*np.sqrt(0+b)

file=open("Raw Data/IMHB/AminoAlcohols/Data.json", 'r')
dataaa=json.loads(file.read())
file.close()
file=open("Raw Data/IMHB/Enolone/Data.json", 'r')
dataen=json.loads(file.read())
file.close()
file=open("Raw Data/Normal/Data.json", 'r')
datanorm=json.loads(file.read())
file.close()
file=open("Raw Data/Methyl/Data.json", 'r')
datameth=json.loads(file.read())
file.close()
file=open("Raw Data/Acids/Data.json", 'r')
dataacid=json.loads(file.read())
file.close()

distshift=0.011
freqshift=150

colorlist=[]

dwhbae=[]
dvhbae=[]
drhbae=[]
dwmae=[]
drmae=[]
deae=[]
deohae=[]
colorsi=[]
for key in dataaa.keys():
    don='h2o'
    basis='avdz'
    if key=="9": continue
    deae.append(dataaa[key]['dEe'])
    deohae.append(dataaa[key]['dEoh'])
    drhbae.append(dataaa[key]['dre'])
    dwhbae.append(dataaa[key]['dw'])
    drmae.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataaa[key]['rm'])
    dwmae.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataaa[key]['wm'])
    colorsi.append('tab:blue')

dwhbee=[]
dvhbee=[]
drhbee=[]
dwmee=[]
dvmee=[]
drmee=[]
deee=[]
deohee=[]
colorsee=[]
for key in dataen.keys():
    don='h2o'
    basis='avdz'
    if key=="9": continue
    deee.append(dataen[key]['dEe'])
    deohee.append(dataen[key]['dEoh'])
    drhbee.append(dataen[key]['dre'])
    dwhbee.append(dataen[key]['dw'])
    drmee.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataen[key]['rm']-distshift)
    dwmee.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataen[key]['wm']+freqshift)
    colorsi.append('tab:red')

dwhbi=np.concatenate((dwhbae,dwhbee))
drhbi=np.concatenate((drhbae,drhbee))
dwmi=np.concatenate((dwmae,dwmee))
drmi=np.concatenate((drmae,drmee))
dei=np.concatenate((deae,deee))
deohi=np.concatenate((deohae,deohee))

dwhb=[]
dvhb=[]
drhb=[]
dwm=[]
dvm=[]
drm=[]
de=[]
deoh=[]
deoa=[]
colors=[]
for key in datanorm.keys():
    de.append(datanorm[key]['dEe'])
    deoh.append(datanorm[key]['dEoh'])
    deoa.append(datanorm[key]['dEoa'])
    drhb.append(datanorm[key]['dre'])
    dwhb.append(datanorm[key]['dw'])
    dvhb.append(datanorm[key]['dv'])
    drm.append(0)
    dwm.append(0)
    dvm.append(0)
    colors.append('tab:blue')
for key in datameth.keys():
    if datameth[key]['donor']=='az' and datameth[key]['acceptor']=='h2o': continue
    if datameth[key]['donor']=='meoh': don='h2o'
    else: don='nh3'
    basis=datameth[key]['basis']
    de.append(datameth[key]['dEe'])
    deoh.append(datameth[key]['dEoh'])
    deoa.append(datameth[key]['dEoa'])
    drhb.append(datameth[key]['dre'])
    dwhb.append(datameth[key]['dw'])
    dvhb.append(datameth[key]['dv'])
    drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-datameth[key]['donor_r'])
    dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-datameth[key]['donor_w'])
    dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-datameth[key]['donor_v'])
    colors.append('tab:brown')
for key in dataacid.keys():
    don='h2o'
    basis=dataacid[key]['basis']
    de.append(dataacid[key]['dEe'])
    deoh.append(dataacid[key]['dEoh'])
    deoa.append(dataacid[key]['dEoa'])
    drhb.append(dataacid[key]['dre'])
    dwhb.append(dataacid[key]['dw'])
    dvhb.append(dataacid[key]['dv'])
    drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataacid[key]['donor_r'])
    dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataacid[key]['donor_w'])
    dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-dataacid[key]['donor_v'])
    colors.append('tab:red')
file=open("Raw Data/ParameterizeGroups/enol/Data.json", 'r')
dataenol=json.loads(file.read())
file.close()
for key in dataenol.keys():
    don='h2o'
    basis=dataenol[key]['basis']
    de.append(dataenol[key]['dEe'])
    deoh.append(dataenol[key]['dEoh'])
    deoa.append(dataenol[key]['dEoa'])
    drhb.append(dataenol[key]['dre'])
    dwhb.append(dataenol[key]['dw'])
    dvhb.append(dataenol[key]['dv'])
    drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataenol[key]['donor_r']-distshift)
    dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataenol[key]['donor_w']+freqshift)
    dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-dataenol[key]['donor_v'])
    colors.append('tab:green')
file = open("Raw Data/ParameterizeGroups/enolone/Data.json", 'r')
dataenolone = json.loads(file.read())
file.close()
for key in dataenolone.keys():
    don='h2o'
    basis=dataenolone[key]['basis']
    de.append(dataenolone[key]['dEe'])
    deoh.append(dataenolone[key]['dEoh'])
    deoa.append(dataenolone[key]['dEoa'])
    drhb.append(dataenolone[key]['dre'])
    dwhb.append(dataenolone[key]['dw'])
    dvhb.append(dataenolone[key]['dv'])
    drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataenolone[key]['donor_r']-distshift)
    dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataenolone[key]['donor_w']+freqshift)
    dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-dataenolone[key]['donor_v']+freqshift)
    colors.append('tab:green')
file=open("Raw Data/ParameterizeGroups/carboxylicAcid/Data.json", 'r')
datacooh=json.loads(file.read())
file.close()
for key in datacooh.keys():
    don='h2o'
    if datacooh[key]['donor']=='acohc' or datacooh[key]['donor']=='hcoohc': continue
    basis=datacooh[key]['basis']
    de.append(datacooh[key]['dEe'])
    deoh.append(datacooh[key]['dEoh'])
    deoa.append(datacooh[key]['dEoa'])
    drhb.append(datacooh[key]['dre'])
    dwhb.append(datacooh[key]['dw'])
    dvhb.append(datacooh[key]['dv'])
    drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-datacooh[key]['donor_r']-distshift)
    dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-datacooh[key]['donor_w']+freqshift)
    dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-datacooh[key]['donor_v']+freqshift)
    colors.append('tab:orange')
file=open("Raw Data/ParameterizeGroups/amide/Data.json", 'r')
dataamide=json.loads(file.read())
file.close()
for key in dataamide.keys():
    don='nh3'
    if dataamide[key]['donor']=='acohc' or dataamide[key]['donor']=='hcoohc': continue
    basis=dataamide[key]['basis']
    de.append(dataamide[key]['dEe'])
    deoh.append(dataamide[key]['dEoh'])
    deoa.append(dataamide[key]['dEoa'])
    drhb.append(dataamide[key]['dre'])
    dwhb.append(dataamide[key]['dw'])
    dvhb.append(dataamide[key]['dv'])
    drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataamide[key]['donor_r']-distshift)
    dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataamide[key]['donor_w']+freqshift)
    dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-dataamide[key]['donor_v']+freqshift)
    colors.append('tab:pink')

file=open("Raw Data/ParameterizeGroups/amide/Data.json", 'r')
dataamide=json.loads(file.read())
file.close()
file=open("Raw Data/ParameterizeGroups/imine/Data.json", 'r')
dataimine=json.loads(file.read())
file.close()
for key in dataamide.keys():
    don='nh3'
    if dataamide[key]['donor']=='aconh2c' or dataamide[key]['donor']=='hconh2c': continue
    basis=dataamide[key]['basis']
    de.append(dataamide[key]['dEe'])
    deoh.append(dataamide[key]['dEoh'])
    deoa.append(dataamide[key]['dEoa'])
    drhb.append(dataamide[key]['dre'])
    dwhb.append(dataamide[key]['dw'])
    dvhb.append(dataamide[key]['dv'])
    drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataamide[key]['donor_r']-distshift)
    dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataamide[key]['donor_w']+freqshift)
    dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-dataamide[key]['donor_v']+freqshift)
    colors.append('tab:pink')

for key in dataimine.keys():
    don='nh3'
    if dataimine[key]['donor']=='aconh2c' or dataimine[key]['donor']=='hconh2c': continue
    basis=dataimine[key]['basis']
    de.append(dataimine[key]['dEe'])
    deoh.append(dataimine[key]['dEoh'])
    deoa.append(dataimine[key]['dEoa'])
    drhb.append(dataimine[key]['dre'])
    dwhb.append(dataimine[key]['dw'])
    dvhb.append(dataimine[key]['dv'])
    if dataimine[key]['donor']=='imine':
        drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataimine[key]['donor_r'])
        dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataimine[key]['donor_w'])
        dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-dataimine[key]['donor_v'])
    else:
        drm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_r']-dataimine[key]['donor_r']-distshift)
        dwm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_w']-dataimine[key]['donor_w']+freqshift)
        dvm.append(datanorm["{}_{}_{}".format(don,don,basis)]['donor_v']-dataimine[key]['donor_v']+freqshift)
    colors.append('tab:purple')

dw=np.array(dwhb)
dr=np.array(drhb)
de=np.array(de)
deoh=np.array(deoh)
deoa=np.array(deoa)
dwm=np.array(dwm)
drm=np.array(drm)*-100
dvm=np.array(dvm)

dwi=np.array(dwhbi)
dri=np.array(drhbi)
dei=np.array(dei)
deohi=np.array(deohi)
dwmi=np.array(dwmi)
drmi=np.array(drmi)*-100

dwae=np.array(dwhbae)
drae=np.array(drhbae)
deae=np.array(deae)
deohae=np.array(deohae)
dwmae=np.array(dwmae)
drmae=np.array(drmae)*-100

dwee=np.array(dwhbee)
dree=np.array(drhbee)
deee=np.array(deee)
deohee=np.array(deohee)
dwmee=np.array(dwmee)
drmee=np.array(drmee)*-100

rs=np.linspace(0,5,10000)
ws=np.linspace(0,1000,10000)
scale=0.025
s=1e-12
plt.rcParams.update({'font.size':10,"font.family":'arial'})
fig,axes=plt.subplots(2,2,sharex='col',sharey='row',figsize=(3,3))
fitre=lm.LinearRegression().fit(np.array([dr,drm]).T,de)
fitreI=opt.minimize(RMSDf,x0=np.array([fitre.coef_[0],0,0.5]),args=(dr,drm,de,Iogansen),
                    bounds=((0,np.inf),(0,np.inf),(0.5-s,0.5+s))).x
print("Ee vs. dr:", fitre.coef_[0],fitre.coef_[1]/fitre.coef_[0],fitre.intercept_)
axes[0,0].plot(rs,rs*fitre.coef_[0]+fitre.intercept_,c='k',linestyle='--',zorder=100)
axes[0,0].plot(rs,Iogansen(rs,0,*fitreI),c='k',zorder=100)
axes[0,0].scatter(dr+0.5*drm,de,c='tab:gray',marker='.',zorder=1,alpha=0.25)
axes[0,0].set_ylim(0,14.95)
axes[0,0].set_xlim(0,4.95)
axes[0,0].scatter(dri+0.5*drmi,dei,c=colorsi,marker='.',zorder=10)
axes[0,0].text(4.95*(scale),14.95*(1-scale),"MSE: {:.2f}\nkcal/mol".format(
    MSE(deae,Iogansen(drhbae,drmae,*fitreI))),
               ha='left',va='top',fontsize=8,c='tab:blue',weight='semibold')
axes[0,0].text(4.95*(1-scale),14.95*(scale),"MSE: {:.2f}\nkcal/mol".format(
    MSE(deee,Iogansen(drhbee,drmee,*fitreI))),
               ha='right',va='bottom',fontsize=8,c='tab:red',weight='semibold')

fitwe=lm.LinearRegression().fit(np.array([dw,dwm]).T,de)
fitweI=opt.minimize(RMSDf,x0=np.array([fitwe.coef_[0],0,1]),args=(dw,dwm,de,Iogansen),
                    bounds=((0,np.inf),(0,np.inf),(1-s,1+s))).x
print("Ee vs. dw:", fitwe.coef_[0],fitwe.coef_[1]/fitwe.coef_[0],fitwe.intercept_)
axes[0,1].set_xlim(0,1000)
axes[0,1].scatter(dw+1*dwm,de,c='tab:gray',marker='.',zorder=1,alpha=0.25)
axes[0,1].scatter(dwi+dwmi,dei,c=colorsi,marker='.',zorder=10)
axes[0,1].plot(ws,ws*fitwe.coef_[0]+fitwe.intercept_,c='k',linestyle='--',zorder=100)
axes[0,1].plot(ws,Iogansen(ws,0,*fitweI),c='k',zorder=100)
axes[0,1].text(1000*(scale),14.95*(1-scale),"MSE: {:.2f}\nkcal/mol".format(
    MSE(deae,Iogansen(dwhbae,dwmae,*fitweI))),
               ha='left',va='top',fontsize=8,c='tab:blue',weight='semibold')
axes[0,1].text(1000*(1-scale),14.95*(scale),"MSE: {:.2f}\nkcal/mol".format(
    MSE(deee,Iogansen(dwhbee,dwmee,*fitweI))),
               ha='right',va='bottom',fontsize=8,c='tab:red',weight='semibold')

fitro=lm.LinearRegression().fit(np.array([dr,drm]).T,deoh)
fitroI=opt.minimize(RMSDf,x0=np.array([fitro.coef_[0],0,0.5]),args=(dr,drm,deoh,Iogansen),
                    bounds=((0,np.inf),(0,np.inf),(0.5-s,0.5+s))).x
print("Eo vs. dr:", fitro.coef_[0],fitro.coef_[1]/fitro.coef_[0],fitro.intercept_)
axes[1,0].set_ylim(0,14.95)
axes[1,0].scatter(dr+0.5*drm,deoh,c='tab:gray',marker='.',alpha=0.25)
axes[1,0].scatter(dri+0.5*drmi,deohi,c=colorsi,marker='.',zorder=10)
axes[1,0].plot(rs,rs*fitro.coef_[0]+fitro.intercept_,c='k',linestyle='--',zorder=100)
axes[1,0].plot(rs,Iogansen(rs,0,*fitroI),c='k',zorder=100)
axes[1,0].text(4.95*(scale),14.95*(1-scale),"MSE: {:.2f}\nkcal/mol".format(
    MSE(deohae,Iogansen(drhbae,drmae,*fitroI))),
               ha='left',va='top',fontsize=8,c='tab:blue',weight='semibold')
axes[1,0].text(4.95*(1-scale),14.95*(scale),"MSE: {:.2f}\nkcal/mol".format(
    MSE(deohee,Iogansen(drhbee,drmee,*fitroI))),
               ha='right',va='bottom',fontsize=8,c='tab:red',weight='semibold')
print(deohee-Iogansen(drhbee,drmee,*fitroI))
fitwo=lm.LinearRegression().fit(np.array([dw,dwm]).T,deoh)
fitwoI=opt.minimize(RMSDf,x0=np.array([fitwo.coef_[0],0,1]),args=(dw,dwm,deoh,Iogansen),
                    bounds=((0,np.inf),(0,np.inf),(1-s,1+s))).x
print("Eo vs. dw:", fitwo.coef_[0],fitwo.coef_[1]/fitwo.coef_[0],fitwo.intercept_)
axes[1,1].set_xlim(0,1000)
axes[1,1].scatter(dw+1*dwm,deoh,c='tab:gray',marker='.',alpha=0.25)
axes[1,1].scatter(dwi+dwmi,deohi,c=colorsi,marker='.',zorder=10)
axes[1,1].plot(ws,ws*fitwo.coef_[0]+fitwo.intercept_,c='k',linestyle='--',zorder=100)
axes[1,1].set_xlabel("$\Delta \omega_\mathrm{HB}+\Delta \omega_\mathrm{M}$\n$(\mathrm{cm}^{-1})$")
axes[1,0].set_xlabel("$\Delta r_\mathrm{HB}+\\frac{1}{2}\Delta r_\mathrm{M}$\n$(\mathrm{pm})$")
axes[0,0].set_ylabel("$\Delta E_e$ (kcal/mol)")
axes[1,0].set_ylabel("$\Delta E_0$ (kcal/mol)")
axes[1,1].plot(ws,Iogansen(ws,0,*fitwoI),c='k',zorder=100)
axes[1,1].text(1000*(scale),14.95*(1-scale),"MSE: {:.2f}\nkcal/mol".format(
    MSE(deohae,Iogansen(dwhbae,dwmae,*fitwoI))),
               ha='left',va='top',fontsize=8,c='tab:blue',weight='semibold')
axes[1,1].text(1000*(1-scale),14.95*(scale),"MSE: {:.2f}\nkcal/mol".format(
    MSE(deohee,Iogansen(dwhbee,dwmee,*fitwoI))),
               ha='right',va='bottom',fontsize=8,c='tab:red',weight='semibold')
plt.tight_layout()
axes[0,1].tick_params(labelleft=False,left=False)
axes[1,1].tick_params(labelleft=False,left=False)
plt.subplots_adjust(hspace=0,wspace=0)
plt.savefig("Figures/IMHBsProve.png",dpi=600)
plt.show()