import numpy as np
import matplotlib.pyplot as plt
import sklearn.linear_model as lm
import scipy.optimize as opt
import json
plt.rcParams.update({'font.size':12,"font.family":'arial'})
fig,axes=plt.subplots(2,2,sharex='col',sharey='row',figsize=(5,5))
def RMSD(e,f):
    return np.sqrt(np.sum((e-f)**2)/len(e))
def RMSDf(params,hb,m,e,func):
    return np.sqrt(np.sum((e-func(hb,m,*params))**2)/len(e))
def Iogansen(hb,m,a,b,c):
    return a*np.sqrt(hb+c*m+b)-a*np.sqrt(0+b)

file=open("Raw Data/Normal/Data.json",'r')
data=json.loads(file.read())
file.close()
file=open("Raw Data/Methyl/Data.json",'r')
datameth=json.loads(file.read())
file.close()
file=open("Raw Data/Acids/Data.json",'r')
dataacid=json.loads(file.read())
file.close()
file=open("Raw Data/ParameterizeGroups/enol/Data.json",'r')
dataenol=json.loads(file.read())
file.close()
file=open("Raw Data/ParameterizeGroups/enolone/Data.json",'r')
dataenolone=json.loads(file.read())
file.close()
file=open("Raw Data/ParameterizeGroups/carboxylicAcid/Data.json",'r')
datacooh=json.loads(file.read())
file.close()
file=open("Raw Data/ParameterizeGroups/amide/Data.json",'r')
dataamide=json.loads(file.read())
file.close()
#"""

file=open("Raw Data/ParameterizeGroups/imine/Data.json",'r')
dataimine=json.loads(file.read())
file.close()#"""

freqshift=150
distshift=0.011

colorlist=[]
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
for key in data.keys():
    de.append(data[key]['dEe'])
    deoh.append(data[key]['dEoh'])
    deoa.append(data[key]['dEoa'])
    drhb.append(data[key]['dre'])
    dwhb.append(data[key]['dw'])
    dvhb.append(data[key]['dv'])
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
    drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-datameth[key]['donor_r'])
    dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-datameth[key]['donor_w'])
    dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-datameth[key]['donor_v'])
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
    drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-dataacid[key]['donor_r'])
    dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-dataacid[key]['donor_w'])
    dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-dataacid[key]['donor_v'])
    colors.append('tab:red')
for key in dataenol.keys():
    don='h2o'
    basis=dataenol[key]['basis']
    de.append(dataenol[key]['dEe'])
    deoh.append(dataenol[key]['dEoh'])
    deoa.append(dataenol[key]['dEoa'])
    drhb.append(dataenol[key]['dre'])
    dwhb.append(dataenol[key]['dw'])
    dvhb.append(dataenol[key]['dv'])
    drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-dataenol[key]['donor_r']-distshift)
    dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-dataenol[key]['donor_w']+freqshift)
    dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-dataenol[key]['donor_v']+freqshift)
    colors.append('tab:green')
for key in dataenolone.keys():
    don='h2o'
    basis=dataenolone[key]['basis']
    de.append(dataenolone[key]['dEe'])
    deoh.append(dataenolone[key]['dEoh'])
    deoa.append(dataenolone[key]['dEoa'])
    drhb.append(dataenolone[key]['dre'])
    dwhb.append(dataenolone[key]['dw'])
    dvhb.append(dataenolone[key]['dv'])
    drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-dataenolone[key]['donor_r']-distshift)
    dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-dataenolone[key]['donor_w']+freqshift)
    dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-dataenolone[key]['donor_v']+freqshift)
    colors.append('tab:green')


for key in datacooh.keys():
    don='h2o'
    if datacooh[key]['donor']=='acohc' or datacooh[key]['donor']=='hcoohc':
        continue
        drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-datacooh[key]['donor_r'])
        dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-datacooh[key]['donor_w'])
        dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-datacooh[key]['donor_v'])
    else:
        drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-datacooh[key]['donor_r']-distshift)
        dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-datacooh[key]['donor_w']+freqshift)
        dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-datacooh[key]['donor_v']+freqshift)
    basis=datacooh[key]['basis']
    de.append(datacooh[key]['dEe'])
    deoh.append(datacooh[key]['dEoh'])
    deoa.append(datacooh[key]['dEoa'])
    drhb.append(datacooh[key]['dre'])
    dwhb.append(datacooh[key]['dw'])
    dvhb.append(datacooh[key]['dv'])
    colors.append('tab:orange')
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
    drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-dataamide[key]['donor_r']-distshift)
    dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-dataamide[key]['donor_w']+freqshift)
    dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-dataamide[key]['donor_v']+freqshift)
    colors.append('tab:pink')
#"""
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
        drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-dataimine[key]['donor_r'])
        dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-dataimine[key]['donor_w'])
        dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-dataimine[key]['donor_v'])
    else:
        drm.append(data["{}_{}_{}".format(don,don,basis)]['donor_r']-dataimine[key]['donor_r']-distshift)
        dwm.append(data["{}_{}_{}".format(don,don,basis)]['donor_w']-dataimine[key]['donor_w']+freqshift)
        dvm.append(data["{}_{}_{}".format(don,don,basis)]['donor_v']-dataimine[key]['donor_v']+freqshift)
    colors.append('tab:purple')
#"""
dw=np.array(dwhb)
dr=np.array(drhb)
de=np.array(de)
deoh=np.array(deoh)
deoa=np.array(deoa)
dwm=np.array(dwm)
drm=np.array(drm)*-100
dvm=np.array(dvm)

"""dwhbe=[]
dvhbe=[]
drhbe=[]
dwme=[]
dvme=[]
drme=[]
dee=[]
deohe=[]
deoae=[]
dwhbc=[]
dvhbc=[]
drhbc=[]
dwmc=[]
dvmc=[]
drmc=[]
dec=[]
deohc=[]
deoac=[]
drhbe=np.array(drhbe)
dwhbe=np.array(dwhbe)
dvhbe=np.array(dvhbe)
drme=np.array(drme)*-100
dwme=np.array(dwme)
dvme=np.array(dvme)
dee=np.array(dee)
drhbc=np.array(drhbc)
dwhbc=np.array(dwhbc)
dvhbc=np.array(dvhbc)
drmc=np.array(drmc)*-100
dwmc=np.array(dwmc)
dvc=np.array(dvmc)

dw=np.concatenate((dw,dwhbe))
dr=np.concatenate((dr,drhbe))
de=np.concatenate((de,dee))
deoh=np.concatenate((deoh,deohe))
deoa=np.concatenate((deoa,deoae))
dwm=np.concatenate((dwm,dwme))
drm=np.concatenate((drm,drme))
dvm=np.concatenate((dvm,dvme))"""

def RMSD2(params,hb,m,e1,e2,func):
    return np.sqrt((np.sum((e1-func(hb,m,params[0],params[2],params[4]))**2)
                    +np.sum((e2-func(hb,m,params[1],params[3],params[4]))**2))/(len(e1)+len(e2)))
def Iogansen(hb,m,a,b,c):
    return a*np.sqrt(hb+c*m+b)-a*np.sqrt(0+b)
s=1e-12
fitrI=opt.basinhopping(RMSD2,x0=np.array([7.0,6.4,0.007,0.09,0.29]),niter=100,
                        minimizer_kwargs={'args':(dr,drm,de,deoh,Iogansen),
                                          'bounds':((0,np.inf),(0,np.inf),(0,np.inf),(0,np.inf),(0.5-s,0.5+s),)}
                        ).x
fitwI=opt.basinhopping(RMSD2,x0=np.array([0.49,0.46,5.3,36,.681]),niter=100,
                        minimizer_kwargs={'args':(dw,dwm,de,deoh,Iogansen),
                                          'bounds':((0,np.inf),(0,np.inf),(0,np.inf),(0,np.inf),(1-s,1+s),)}
                        ).x

rs=np.linspace(0,5,1000)
ws=np.linspace(0,1000,1000)
scale=0.025
#np.set_printoptions(suppress=True,precision=4)
fitre=lm.LinearRegression().fit(np.array([dr+0.5*drm]).T,de)
fitreI=np.array([fitrI[0],fitrI[2],fitrI[4]])
print("Ee vs. dr:", fitre.coef_[0],fitre.intercept_)
print(fitreI)
axes[0,0].set_ylim(0,14.5)
axes[0,0].set_xlim(0,4.5)
axes[0,0].scatter(dr+fitreI[2]*drm,de,c=colors,marker='.')
axes[0,0].plot(rs,rs*fitre.coef_[0]+fitre.intercept_,c='k',linestyle='--')
axes[0,0].plot(rs,Iogansen(rs,0,*fitreI),c='k')
axes[0,0].text(4.5*(1-scale),14.5*scale,"RMSD:\n{:.2f} kcal/mol".format(RMSD(de,Iogansen(dr,drm,*fitreI))),
               ha='right',va='bottom',fontsize=10)

fitwe=lm.LinearRegression().fit(np.array([dw+dwm]).T,de)
fitweI=np.array([fitwI[0],fitwI[2],fitwI[4]])
print("Ee vs. dw:", fitwe.coef_[0],fitwe.intercept_)
print(fitweI)
axes[0,1].set_xlim(0,1000)
axes[0,1].scatter(dw+fitweI[2]*dwm,de,c=colors,marker='.')
axes[0,1].plot(ws,ws*fitwe.coef_[0]+fitwe.intercept_,c='k',linestyle='--')
axes[0,1].plot(ws,Iogansen(ws,0,*fitweI),c='k')
axes[0,1].text(1000*(1-scale),14.5*scale,"RMSD:\n{:.2f} kcal/mol".format(RMSD(de,Iogansen(dw,dwm,*fitweI))),
               ha='right',va='bottom',fontsize=10)

fitro=lm.LinearRegression().fit(np.array([dr+0.5*drm]).T,deoh)
fitroI=np.array([fitrI[1],fitrI[3],fitrI[4]])
print("Eo vs. dr:", fitro.coef_[0],fitro.intercept_)
print(fitroI)
axes[1,0].set_ylim(0,14.5)
axes[1,0].scatter(dr+fitroI[2]*drm,deoh,c=colors,marker='.')
print(len(deoh))
axes[1,0].plot(rs,rs*fitro.coef_[0]+fitro.intercept_,c='k',linestyle='--')
axes[1,0].plot(rs,Iogansen(rs,0,*fitroI),c='k')
axes[1,0].text(4.5*(1-scale),14.5*scale,"RMSD:\n{:.2f} kcal/mol".format(RMSD(deoh,Iogansen(dr,drm,*fitroI))),
               ha='right',va='bottom',fontsize=10)

fitwo=lm.LinearRegression().fit(np.array([dw+dwm]).T,deoh)
fitwoI=np.array([fitwI[1],fitwI[3],fitwI[4]])
print("Eo vs. dw:", fitwo.coef_[0],fitwo.intercept_)
print(fitwoI)
axes[1,1].set_xlim(0,1000)
axes[1,1].scatter(dw+fitwoI[2]*dwm,deoh,c=colors,marker='.')
axes[1,1].plot(ws,ws*fitwo.coef_[0]+fitwo.intercept_,c='k',linestyle='--')
axes[1,1].plot(ws,Iogansen(ws,0,*fitwoI),c='k')
axes[1,1].text(1000*(1-scale),14.5*scale,"RMSD:\n{:.2f} kcal/mol".format(RMSD(deoh,Iogansen(dw,dwm,*fitwoI))),
               ha='right',va='bottom',fontsize=10)

#print(acidic)
#acidic,_=opt.curve_fit(lambda x, a,b: b*x+a, drma,dema)
#print(acidic)


"""axes[0,0].scatter(drhbe+fitre.coef_[1]/fitre.coef_[0]*drme,dee,c='tab:green',marker='.')
axes[0,1].scatter(dwhbe+fitwe.coef_[1]/fitwe.coef_[0]*dwme,dee,c='tab:green',marker='.')
axes[1,0].scatter(drhbe+fitro.coef_[1]/fitro.coef_[0]*drme,deohe,c='tab:green',marker='.')
axes[1,1].scatter(dwhbe+fitwo.coef_[1]/fitwo.coef_[0]*dwme,deohe,c='tab:green',marker='.')"""

"""axes[0,0].scatter(drhbc+fitre.coef_[1]/fitre.coef_[0]*drmc,dec,c='tab:orange',marker='.')
axes[0,1].scatter(dwhbc+fitwe.coef_[1]/fitwe.coef_[0]*dwmc,dec,c='tab:orange',marker='.')
axes[1,0].scatter(drhbc+fitro.coef_[1]/fitro.coef_[0]*drmc,deohc,c='tab:orange',marker='.')
axes[1,1].scatter(dwhbc+fitwo.coef_[1]/fitwo.coef_[0]*dwmc,deohc,c='tab:orange',marker='.')"""

axes[0,0].set_ylabel("$\\Delta E_e$ (kcal/mol)")
axes[1,0].set_ylabel("$\\Delta E_0$ (kcal/mol)")
axes[1,0].set_xlabel("$\Delta r_\\mathrm{HB}+\\frac{1}{2}\\Delta r_\\mathrm{M}$\n(pm)")

fitro=lm.LinearRegression().fit(np.array([dr,drm]).T,deoh)

axes[1,1].set_xlabel("$\Delta \\omega_\\mathrm{HB}+\\Delta \\omega_\\mathrm{M}$\n(cm$^{-1}$)")
fitwo=lm.LinearRegression().fit(np.array([dw,dwm]).T,deoh)



plt.tight_layout()
plt.subplots_adjust(hspace=0,wspace=0)
plt.savefig("Figures/2DFitIFix.png",dpi=600)
plt.show()

aedw=3706-3570
print("Aminoethanol:"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3706-3530
print("Bimethyl Aminoethanol:"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3706-3333
print("Bi-trisfluoro-methyl Aminoethanol:"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
print()
aedw=3407-3298+freqshift
print("C-term"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3407-3320+freqshift
print("Internal"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3407-3351+freqshift
print("Internal"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3407-3362+freqshift
print("Internal"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3407-3404+freqshift
print("NH-pi"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3407-3434+freqshift
print("N-term"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
print()
aedw=3407-(3363+3389)/2+freqshift
print("Beta"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
print()
aedw=3407-3130+freqshift
print("CG"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3407-3165+freqshift
print("CG"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3407-3290+freqshift
print("CG"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
print()
aedw=3407-2784+freqshift
print("AT"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)
aedw=3407-3205+freqshift
print("AT"
      ,aedw,"cm^-1",
      Iogansen(aedw,0,*fitweI),"kcal/mol",
      Iogansen(aedw,0,*fitwoI),"kcal/mol",)