import numpy as np
import math
import matplotlib.pyplot as plt

class Value:
    def __init__(self, list1, list2):
        self.list1 = list1
        self.list2 = list2
        self.mean_x = sum([x for x in self.list1])/len(self.list1)
        self.mean_y = sum([y for y in self.list2])/len(self.list2)
        self.mean_x2 = sum([x**2 for x in self.list1])/len(self.list1)
        self.mean_y2 = sum([y**2 for y in self.list2])/len(self.list2)
        xyl = []
        for i in range(len(self.list1)):
            xyl.append(self.list1[i]*self.list2[i])
        self.mean_xy = sum([xy for xy in xyl])/len(xyl)

    def mean_v(self):
        a = (self.mean_xy - self.mean_x*self.mean_y)/(self.mean_x2 - (self.mean_x)**2)
        b = self.mean_y - a*self.mean_x
        return a, b
    
    def devijacija_a_b(self): #devijacija za y = ax+b
        a, b = self.mean_v()
        sigma_a = np.sqrt((1/len(self.list1))*((self.mean_y2 - self.mean_y**2)/(self.mean_x2 - self.mean_x**2)- a**2))
        sigma_b = sigma_a*np.sqrt(self.mean_x2 - self.mean_x**2)
        return sigma_a, sigma_b
    
    def mean_v2(self): #y =ax
        a = (self.mean_xy)/(self.mean_x2)
        return a
    
    def mean_a(self): #devijacija za oblik y=ax
        val = self.mean_v2()
        dev = np.sqrt((1/len(self.list1)*((self.mean_y2/self.mean_x2) - val**2)))
        return dev
        

class Aritm:
    def __init__(self, lista):
        self.listt = lista
    def sredina(self):
        return sum([t for t in self.listt])/len(self.listt)
    
class Dev:
    def __init__(self, listt):
        self.listt = listt
    def sredina(self):
        return sum([t for t in self.listt])/len(self.listt)
    def devijacija(self):
        brojnik = 0
        mean_t = sum([t for t in self.listt])/len(self.listt)
        for el in self.listt:
            brojnik += (el - mean_t)**2
        return np.sqrt(brojnik/(len(self.listt)*(len(self.listt)- 1)))
    


