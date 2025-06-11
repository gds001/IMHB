import numpy as np
import matplotlib.pyplot as plt
def RMSD(e,f):
    return np.sqrt(np.sum((e-f)**2)/len(e))
def MSE(e,f):
    return np.sum((e-f))/len(e)
def MAE(e,f):
    return np.sum(np.abs(e-f))/len(e)
plt.rcParams.update({'font.size':12,'font.family':'arial'})
data=np.genfromtxt('Processed Data/Anharms.csv', delimiter=',').T

plt.figure(figsize=(4,4))
plt.scatter(data[0],data[2]-data[0],marker='.')
plt.scatter(data[0],data[1]-data[0],zorder=0,c='tab:red',marker='.')
#plt.legend(["$\\Delta \\nu_\mathrm{HB}$","$\Delta \\nu_\mathrm{Dep.}$"])
plt.plot([0,1000],[0,0],c='k')
plt.xlabel("$\\Delta \omega_\mathrm{HB}$ (cm$^{-1}$)")
plt.ylabel("$\\Delta \\nu_\mathrm{HB}-\Delta \omega_\mathrm{HB}$ (cm$^{-1}$)")
plt.text(1000*0.975,1100*0.975,"RMSD: 32 cm$^{-1}$\nMAD: 27 cm$^{-1}$\nMSD: -13 cm$^{-1}$",ha='right',va='top',color='tab:blue',fontweight='semibold')
plt.text(1000*0.975,-1000*0.975,"RMSD: 206 cm$^{-1}$\nMAD: 92 cm$^{-1}$\nMSD: 29 cm$^{-1}$",ha='right',va='bottom',color='tab:red',fontweight='semibold')
print(RMSD(data[2],data[0]),RMSD(data[1],data[0]))
print(MAE(data[2],data[0]),MAE(data[1],data[0]))
print(MSE(data[2],data[0]),MSE(data[1],data[0]))
plt.xlim(0,1000)
plt.ylim(-1000,1100)
plt.tight_layout()
plt.savefig("Figures/Anharms.png",dpi=600)
plt.show()
