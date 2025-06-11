import numpy as np
import matplotlib.pyplot as plt
import sklearn.linear_model as lm
import scipy.optimize as opt
import json
plt.rcParams.update({'font.size':10,"font.family":'arial'})
def RMSD(e,f):
    return np.sqrt(np.sum((e-f)**2)/len(e))
def Iogansen(x,a,b): return a*np.sqrt(x+b)-a*np.sqrt(b)

file=open("Data/Normal/Data.json", 'r')
data=json.loads(file.read())
file.close()

fig,axes=plt.subplots(2,2,sharex='col',sharey='row',figsize=(3,3))

colorlist=[]
dw=[]
dv=[]
dr=[]
de=[]
deoh=[]
deoa=[]
wm=[]
for key in data.keys():
    de.append(data[key]['dEe'])
    deoh.append(data[key]['dEoh'])
    deoa.append(data[key]['dEoa'])
    dr.append(data[key]['dre'])
    dw.append(data[key]['dw'])
    dv.append(data[key]['dv'])
    if data[key]['donor']=='h2o': wm.append(data[key]['donor_w'])
print(np.average(wm))

dw=np.array(dw)
dr=np.array(dr)
de=np.array(de)
deoh=np.array(deoh)
deoa=np.array(deoa)


rs=np.linspace(0,5,1000)
ws=np.linspace(0,1000,1000)
scale=0.025

axes[0,0].scatter(dr,de,marker='.',c='tab:blue')
fitre=lm.LinearRegression().fit(np.array([dr]).T,de)
fitreI,_=opt.curve_fit(Iogansen,dr,de,p0=np.array([fitre.coef_[0],0]))
#acidic,_=opt.curve_fit(lambda x, a: fit.coef_[0]*x+a, drma,dema)
axes[0,0].plot(rs,Iogansen(rs,*fitreI),c='k')
axes[0,0].plot(rs,fitre.coef_[0]*rs+fitre.intercept_,c='k',linestyle='--')
axes[0,0].set_ylim(0,14.5)
axes[0,0].set_xlim(0,4.5)
axes[0,0].set_ylabel("$\\Delta E_e$ (kcal/mol)")


print("Ee vs. dr",fitre.coef_[0],fitre.intercept_)
print(fitreI)
#print(acidic)
#acidic,_=opt.curve_fit(lambda x, a,b: b*x+a, drma,dema)
#print(acidic)

axes[0,1].scatter(dw,de,marker='.',c='tab:blue')
fitwe=lm.LinearRegression().fit(np.array([dw]).T,de)
fitweI,_=opt.curve_fit(Iogansen,dw,de,p0=np.array([fitwe.coef_[0],0]))
#acidic,_=opt.curve_fit(lambda x, a: fit.coef_[0]*x+a, dwma,dema)
axes[0,1].plot(ws,Iogansen(ws,*fitweI),c='k')
axes[0,1].plot(ws,fitwe.coef_[0]*ws+fitwe.intercept_,c='k',linestyle='--')
#axes[0,1].plot(ws,fit.coef_[0]*ws+acidic[0],c='tab:red')
axes[0,1].set_xlim(0,1000)
print("Ee vs. dw",fitwe.coef_[0],fitwe.intercept_)
print(fitweI)
#print(acidic)
#acidic,_=opt.curve_fit(lambda x, a,b: b*x+a, dwma,dema)
#print(acidic)
axes[1,0].set_ylabel("$\\Delta E_0$ (kcal/mol)")
axes[1,0].set_xlabel("$\\Delta r_\\mathrm{HB}$ (pm)")

axes[1,0].scatter(dr,deoh,marker='.',c='tab:blue')
fitro=lm.LinearRegression().fit(np.array([dr]).T,deoh)
#acidic,_=opt.curve_fit(lambda x, a: fit.coef_[0]*x+a, drma,deoma)
fitroI,_=opt.curve_fit(Iogansen,dr,deoh,p0=np.array([fitro.coef_[0],0]))
axes[1,0].plot(rs,Iogansen(rs,*fitroI),c='k')
axes[1,0].plot(rs,fitro.coef_[0]*rs+fitro.intercept_,c='k',linestyle='--')
#axes[1,0].plot(rs,fit.coef_[0]*rs+acidic[0],c='tab:red')
axes[1,0].set_ylim(0,14.5)
axes[1,0].set_xticks([0,2,4],[0,2,4])
print("Eo vs. dr",fitro.coef_[0],fitro.intercept_)
print(fitroI)
#print(acidic)
#acidic,_=opt.curve_fit(lambda x,a,b: b*x+a, drma,deoma)
#print(acidic)

axes[1,1].set_xlabel("$\\Delta \\omega_\\mathrm{HB}$ (cm$^{-1}$)")
axes[1,1].scatter(dw,deoh,marker='.',c='tab:blue')
fitwo=lm.LinearRegression().fit(np.array([dw]).T,deoh)
#acidic,_=opt.curve_fit(lambda x, a: fit.coef_[0]*x+a, dwma,deoma)
fitwoI,_=opt.curve_fit(Iogansen,dw,deoh,p0=np.array([fitwo.coef_[0],0]))
axes[1,1].plot(ws,Iogansen(ws,*fitwoI),c='k')
axes[1,1].plot(ws,fitwo.coef_[0]*ws+fitwo.intercept_,c='k',linestyle='--')
#axes[1,1].plot(ws,fit.coef_[0]*ws+acidic[0],c='tab:red')
print("Eo vs. dw",fitwo.coef_[0],fitwo.intercept_)
print(fitwoI)
#print(acidic)
#acidic,_=opt.curve_fit(lambda x, a,b: b*x+a, dwma,deoma)
#print(acidic)
s=1e-8
plt.tight_layout()
plt.subplots_adjust(hspace=0,wspace=0)
#plt.savefig("Monomer Reference 1.png",dpi=600)
zzp=Iogansen(dr,*fitreI)
zop=Iogansen(dw,*fitweI)
ozp=Iogansen(dr,*fitroI)
oop=Iogansen(dw,*fitwoI)
zzpl=fitre.coef_[0]*dr+fitre.intercept_
zopl=fitwe.coef_[0]*dw+fitwe.intercept_
ozpl=fitro.coef_[0]*dr+fitro.intercept_
oopl=fitwo.coef_[0]*dw+fitwo.intercept_
#plt.savefig("Monomer Reference 2.png",dpi=600)
zze=de
zoe=de
oze=deoh
ooe=deoh
axes[0,0].text(4.5*(1-scale),14.5*scale,"RMSD:\n{:.2f} kcal/mol".format(RMSD(zze,zzp)),
               ha='right',va='bottom',fontsize=9)
axes[0,1].text(1000*(1-scale),14.5*scale,"RMSD:\n{:.2f} kcal/mol".format(RMSD(zoe,zop)),
               ha='right',va='bottom',fontsize=9)
axes[1,0].text(4.5*(1-scale),14.5*scale,"RMSD:\n{:.2f} kcal/mol".format(RMSD(oze,ozp)),
               ha='right',va='bottom',fontsize=9)
axes[1,1].text(1000*(1-scale),14.5*scale,"RMSD:\n{:.2f} kcal/mol".format(RMSD(ooe,oop)),
               ha='right',va='bottom',fontsize=9)
#axes[1,1].plot(ws,0.0124*ws+0.51,c='k')
plt.savefig("Monomer Reference 1 I.png",dpi=600)
plt.show()

fig,axes=plt.subplots(1,1,figsize=(2,2))
axes.set_xlabel("$\\Delta \\omega_\\mathrm{HB}$ (cm$^{-1}$)")
axes.set_ylabel("$\\Delta r\\mathrm{HB}$ (pm)")
axes.scatter(dw,dr,marker='.',c='tab:blue')
fit=lm.LinearRegression(fit_intercept=False).fit(np.array([dw]).T,dr)
dwp=fit.coef_[0]*dw+fit.intercept_
axes.plot(ws,fit.coef_[0]*ws+fit.intercept_,c='k')
print("dr vs. dw",fit.coef_[0],fit.intercept_)

drs=dr

axes.text(1000*(1-scale),4.5*scale,"RMSD:\n{:.2f} pm".format(RMSD(drs,dwp)),
               ha='right',va='bottom',fontsize=10)
axes.plot(ws,fit.coef_[0]*ws+fit.intercept_,c='tab:red')
print("dr vs. dw",fit.coef_[0],fit.intercept_)
plt.xlim(0,1000)
plt.ylim(0,4.25)
plt.tight_layout()
plt.savefig("StructureSpectra 1.png",dpi=600)
plt.close()