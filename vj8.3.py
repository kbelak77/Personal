import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#definiramo vrij u t=0

L= 20 #domena prostora
b= 10**(-2)
Nx = 100
Nt = 500
dx = (L)/(Nx + 1)
dt = 0.5


alpha = b * (dt/(dx**2))
if alpha >=(1/2):
    print('Rješenje nije stabilno za ovaj odabir dx.')
#rubni uvjeti su a=b=0

#funkcija je dvodimenz polje, u prvom redku postavljamo vrijednosti x-a za t=0 kako je zadano

#ja u biti postavljam bvrijednosti u funkciju za odredeni t odnosno x kad mmi je poznat iznos funkcije, matrica je u biti pohrana vrijednosti funkcije za dana cvorista
#potrebna je matr za pohranu tih podataka a pristupam njenim elementima u petlji za rekurzivno racunati iznose funkcije za sve x odn t
#svaki redak je "slika" fje u vremenu t
#stupac u prostoru x
#poznajemo vrij za xeve, trazimo vrij u vremenskoj evoluciji koristeci formulu
#definiramo domene za t odn x

x = np.linspace(0, L , Nx) #u Nx cvorova na ovoj domeni racunamo
t = np.linspace(0, Nt*dt, Nt) #T=Nt*dt je finalno vrijeme
fja = np.zeros((Nt, Nx))
#redak= nulti za poc uvjete odn t=0

for k in range(len(x)): #za t=0 odn pocetne uvjete pridodaje vrijednosti fje koje su zadane, a ovisne o x
    if (x[k] >= 2) and (x[k] <= 5):
        fja[0, k] = 5.5
    else:
        fja[0,k] = 0.0
fja2 = fja.copy()
fja3=fja.copy()
#u matricu fje ne upisujem listu xeva odn t jer to nije lista toga vec poznate vrijednosti fje
#u tim tockama koje poznajemo
#rubni uvjeti
#rekurzivno popunjavanje iznosa funkcije prema formuli uz pomoc poznatih p.u., r.u.
for j in range(0, Nt-1): #radi se o indeksaciji pa krece od nula bit ce viska
    for i in range(1, Nx-1): #krece od x=1 jer je nula za x=0 i za x=20 sto je posljednji clan u x
        fja[j+1, i] = alpha*fja[j, i+1] + (1-2*alpha)*fja[j, i] + alpha*fja[j, i-1]
#iz rubnih uv
#print(fja2)

#Implicitno ovaj mora matrično? ispunjen je uvjet a=b=0
#V je stupčasta matrica sa vrijednostima funkcije u istom vremenu ali varijabilnim prostornim koord
#odnosno Vj je jedan stupac što nosi vrijednosti matrice u vrem trenutku/koraku j
#doznajemo sve vrij za sve t rješavajući rekurz
#obje umanjene za 2 zbog elimin p.u.

A = np.zeros((Nx - 2, Nx - 2))
np.fill_diagonal(A, 1 + 2 * alpha) #popunjavanje matr prema sabloni
np.fill_diagonal(A[1:], -alpha)
np.fill_diagonal(A[:, 1:], -alpha)

#print(A)
#fja2
#indeksacija
#Nepoznanica je Vj+1 B = Vj-1
for ii in range(0, Nt-1):
     #pristupa n-tom redku matrice;vrijeme
     #V je matrica sa elementima iz spec redka za spec t(fiksan t)?
    V_j = fja2[ii, 1:-1] #bez prvog i zadnjeg clana jer su ionako nula za t=0 i t=T
    #Vj-1 = A*Vj da nerac inverz
    fja2[ii-1, 1:-1] = np.linalg.solve(A, V_j) #rac iz prethodnog iduci

#fja3
B = np.zeros((Nx - 2, Nx - 2))
np.fill_diagonal(B, 2) #popunjavanje matr prema sabloni
np.fill_diagonal(B[1:], -1)
np.fill_diagonal(B[:, 1:], -1)

for iii in range(0, Nt-1):
    V_jp = np.dot((2*np.eye(Nx-2) - alpha*B), fja3[iii-1, 1:-1]) #prosli
    fja3[iii, 1:-1] = np.linalg.solve(A,V_jp)


js = [x for x in range(0,500, 100)]
#print(js)
#traži graf fje i x/dx
plt.figure(figsize=(12, 6))
for el in js:
    plt.plot(x/dx, fja[el], label=f"Explic. t={el * dt:.1f}s")
    plt.plot(x/dx, fja2[el], '--', label=f"Implic. t={el * dt:.1f}s")
    plt.plot(x/dx, fja3[el], ':', label=f"CN t={el * dt:.1f}s")

plt.xlabel("x/dx (m)")
plt.ylabel("ρ (kg/m)")
plt.legend()
plt.show()