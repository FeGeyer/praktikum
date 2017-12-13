import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(0, 3, 1000)
me= const.electron_mass
n=3.3543 # Brechungsindex für GaAs bei einer Wellenlänge von 1771.14
A= (const.e**3) * B / ( 8*(np.pi**2) * const.epsilon_0*(const.c**3) *n)  #konst zum bestimmen der eff Masse
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100

# d = dicke= Länge im Strahlengang
hd = 5.11*10**-3
d1 = 1.296*10**-3
d2 = 1.36*10*10**-3
#Werte einelsen l = Wellenläng. tb= teta mit bfeld , to= teta ohne (b/o)bm= bogenmin mit und ohne bfeld
#hochreines GaAg
hl , htb , hbbm , hto , hobm = np.genfromtxt('hochreinesGaAs.txt', unpack = True)
#1. n-dotiertes GaAs
l1 , tb1, bbm1 , to1 , obm1 = np.genfromtxt('1.n-dotiertesGaAs.txt', unpack = True)
#2. n-dotiertes GaAs
l2  , tb2 , bbm2 , to2 , obm2 = np.genfromtxt('2.n-dotiertesGaAs.txt', unpack = True)
#jetz noch die Winkel um und zusammenrechnen. dann die ohne bFeld von den mit BFeld abziehen. Alle Messwerte in ein Plot mit 'o-'
#Winkel zusammenrechnen eine bogenmin = (1°/60) Grad
u = 1/60 #u der kleine Umrechnungsfaktor, weil faul
htb+= hbbm*u
hto+= hobm*u
tb1+= bbm1*u
to1+= obm1*u
tb2+= bbm2*u
to2+= obm2*u
# The True Teta ausrechnen t mit BFeld - t ohne
ht = 0.5*(htb - hto)
t1 = 0.5*(tb1 - to1)
t2 = 0.5*(tb2 - to2)
ht = abs(ht)
t1 = abs(t1)
t2 = abs(t2)
#Umrechnen der Wellenlängen
#hl*=10**-6
#l1*=10**-6
#l2*=10**-6
#Differenzzwischen der hochreinen und den zwei Proben
D1 = abs(ht/hd-t1/d1)
D2 = abs(ht/hd-t2/d2)
#D = np.array([D1[0:1]])
#D = np.append(D,D1[2:])
#L = np.array([hl[0:1]])
#L = np.append(L,hl[2:])
print('Mittelwert der Wellenlängen:' , np.mean(hl))
#Fit
C1= A*2.8*10**16 #weil umgerechnet auf Meter
C2= A*1.2*10**16
def fitf1(x,a,c):
	return x*a + c
def fitf2(x,b, d):
	return x*b + d
D = np.array(D2[0:3])
D = np.append(D,D2[7:])
h = np.array(hl[0:3])
h = np.append(h,hl[7:])
print(D2,hl)
print(D,h)
# params , cov = curve_fit(fitf1 ,hl**2,D1)
# params = correlated_values(params, cov)
# params2 , cov2 = curve_fit(fitf2 ,h,D )
# params2 = correlated_values(params2, cov2)
# a = params[0]
# c = params[1]
# b = params2[0]
# d = params2[1]
# print(a*10**12,c)
# print(b*10**12,d)
# print('Hildegard1:', unp.sqrt(C1/(a*10**12)))
# print('Hildegard2:', unp.sqrt(C2/(b*10**12)))
# print('Hildegard1/me:', unp.sqrt((C1/(a*10**12)))/me)
# print('Hildegard2/me:', unp.sqrt((C2/(b*10**12)))/me)
# H = np.mean([unp.sqrt((C1/(a*10**12)))/me,unp.sqrt((C2/(b*10**12)))/me])
# print('Mittlere Hildegard:' ,H)
# print('Relf:',relf(0.067,H))
# #Tabelle
# #np.savetxt('HGaAstab.txt',np.column_stack([(hl*10**(-6)),ht,(ht/hd)]), delimiter=' & ',newline= r'\\'+'\n' )
# #np.savetxt('1.nGaAstab.txt',np.column_stack([(l1*10**(-6)),t1,(t1/d1)]), delimiter=' & ',newline= r'\\'+'\n' )
# np.savetxt('2.nGaAstab.txt',np.column_stack([(l2*10**(-6)),t2,(t2/d2)]), delimiter=' & ',newline= r'\\'+'\n' )
#
# plt.plot(hl**2 , D1, 'ro', label= 'Differenz mit 2.8.n-dot')
# plt.plot(x**2,fitf1(x**2,noms(a),noms(c)),'r--', label = 'Ausgleichsgerade für 2.8 n-dotiert')
# plt.plot(x**2,fitf2(x**2,noms(b), noms(d)),'b--', label = 'Ausgleichsgerade für 1.2 n-dotiert')
# plt.plot(h**2, D, 'bo', label='Differenz mit 1.2.n-dot')
# plt.plot(hl[3:7]**2, D2[3:7], 'kx', label = 'Nicht betrachtete Werte')
# plt.xlabel(r'Wellenlänge zum Quadrat $\lambda^2 \:/\:\mu m^2 $')
# plt.ylabel(r'Differenz der Längennormierten Winkel $\frac{\Delta\theta}{L} / \degree m^{-1} $')
# plt.legend(loc='best')
# #plt.show()
# plt.savefig('DeltaTheta.pdf')
# #plt.title('Farady-Rotation der GaAs-Proben im Vergleich')
# #plt.plot((l2**2),t2/d2,'go',label='n-dotiertes GaAs mit $N=1,2\cdot10^{18}cm^{-3}FF$')
# #plt.plot((l1**2),t1/d1,'bo',label=r'n-dotiertes GaAs mit $ N=2,8\cdot10^{18}cm^{-3}$')
# #plt.plot((hl**2),ht/hd,'ro', label=r'Hochreines GaAs')
#plt.xlabel(r'Wellenlänge zum Quadrat $\lambda^2 \:/\: m^2 $')
#plt.ylabel(r'Winkelrotation normiert auf Probenlänge $\frac{\theta}{L} \:/\:\degree m^{-1}$')
#plt.legend(loc='best')
#plt.savefig('GaAsimVgl.pdf')
#plt.show()
