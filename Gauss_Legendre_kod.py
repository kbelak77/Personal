import math
import numpy as np

def gauleg(x1, x2, n):
    
 #   Za zadanu donju i gornju granicu integracije x1 i x2, i stupanj polinoma na ova funkcija vraća listu x i w duljine n, koji sadrže nultočke 
 #   i težine za kvadratur Gauss-Legendre kvadraturu na zadanom intervalu.  
  
    EPS = 1e-6  # Preciznost 
    m = (n + 1) // 2 #jedna nultočka za nep n je nula, ostale simetrične--> tražimo pola nult
    #interval def lGR pol -1,1 --> pomocu ovog pozicioniramo nult LP na interval a,b
    #transf interval it -1,1 na a,b
    xm = 0.5 * (x2 + x1) #sredina int
    xl = 0.5 * (x2 - x1) #pola dulj int
    #da tocke na krajevima u int -1,1
    # Inicijalizacija
    #za zapis x,w
    x = np.zeros(n)
    w = np.zeros(n)
    #trazenje nultocaka fje z
    for i in range(1, m + 1): #jedno izvrsavanje jedna nultocka L polinoma
        z = math.cos(math.pi * (i - 0.25) / (n + 0.5)) #prva procjena nultocke, eksperimentalno određena forma
        while True: 
            p1 = 1.0 #za Legr polinom poc uvjet Legr. polinoma u pretpostavljenoj nultocki z mora dat 1
            p2 = 0.0 # nepostojeći--> zbog potreba rekurzije #vrijednost prethodnog polinoma #ALI nije stalno nula prolaskom kroz petlju poprima nenula vrijednosti; bitno samo incijalno da je nula

            for j in range(1, n + 1): #j+1 red l polinoma--> rekurzivno racunamo svaki prisutni red LP
                #rekurzivno racunanje polinoma
                p3 = p2
                p2 = p1
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j #Vrijednost trenutnog Lpn(z)
            pp = n * (z * p1 - p2) / (z * z - 1.0) #rekurz izr derivacija LPn(z)
            z1 = z #pohrana inicijalne pretp
            z = z1 - p1 / pp #nova pretpostavka nultočke NR METODA sad kad imamo L i der L
            if abs(z - z1) < EPS: #uvjet Newt-R metode
                break
        x[i - 1] = xm - xl * z #uz određenu jednu nultočku z za i-to izvrs petlje
        x[n - i] = xm + xl * z #posljednji
        w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp) #tezine racunamo prema formuli
        w[n - i] = w[i - 1] #simetrija
    
    return x, w

