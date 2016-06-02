import numpy as np
import os
import matplotlib.pyplot as plt


i = 0
rho_vect=[.0,.3,.6,.9]
p_vect=[.001,.01,.1]
test_num = 3
cd = os.getcwd()
plt.close()

#rho = 0
#for p in p_vect:
d = [np.array([[0,1],[2,3],[4,1],[6,1]]), np.array([[0,1],[1,2],[3,2],[5,1]])]
data = []
t = []
times = []
Ivt = []
for i in range(2):
#    data = np.loadtxt("%s/Box Sync/Ebola MCM2015/Code/N_50_Infection_v_Time/rho=0_p=0.01_N=50_K=1_test%g.csv" % (cd, i), delimiter = ",")
    data = d[i]
    Ivt += list(np.transpose(data))
    t += list(Ivt[0])
times = list(set(t))
times.sort()
I_avg = range(len(times))
for i in range(len(times)-1): # for each time frame
#    reached_end = False
    for k in range(len(Ivt[0])):
#        if i == len(times)-1:
 #           reached_end = True
        if (times[i] <= Ivt[0][k]) and (Ivt[0][k] < times[i+1]):
            I_avg[i] += Ivt[0][k]
            break
        
times = np.array(times)    
I_avg = np.array(I_avg)
I_avg /= test_num
        

plt.plot(times, I_avg)

plt.xlabel("$t$")
plt.ylabel("$I$")
plt.show()