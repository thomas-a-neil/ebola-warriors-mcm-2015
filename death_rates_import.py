import matplotlib.pyplot as plt
import numpy as np
import os

curDir = os.getcwd()
N = 50
K = 1

rho_vec = np.append(np.append(np.linspace(0, 0.3, 3), np.linspace(0.35, 0.65, 10)), np.linspace(0.7, 0.95, 3))
p_vec = np.logspace(-1, 0, 5)

##to_box = "Box Sync/Ebola MCM2015/Code"
#folder = "N_50_death_rates" #+ "/randomvacc"
#sub_p_v = [0.5]
#
#types = ["g", "r", "b+", "bo", "r"]
#end = len(rho_vec)
#i = 0
#
#sub = [p_vec[0], p_vec[-1]]
#for p in sub:
#    avg_r = []
#    for rho in rho_vec[0:end]:
#        path = os.path.normpath("%s/%s/pr_death_rates_rho=%g_p=%g_N=%g_K=%g.csv" % (curDir, folder, rho, p, N, K))
#        r_v = np.loadtxt(path, delimiter = " ")
#        avg_r += [np.average(r_v)]
#    plt.plot(rho_vec[0:end], avg_r, types[i])
#    i += 1
#
#plt.xlabel("$\\rho$")
#plt.ylabel("$r$")
##plt.title("Death rate $r$ vs. $\\rho$ for $p$ = 0.1 (green) and $p$ = 1 (red)")
#plt.show()
        


# I Want to take p =1, and compare the plots for random and non_random vaccination

sub_p_v = [0.5]

types = ["g", "r", "b+", "bo", "r"]



i = 0
for p in p_vec:
    folder = "N_50_death_rates" + "/random_vacc"
    end = len(rho_vec)
    avg_r = []
    plt.figure(i)
    for rho in rho_vec[0:end]:
        path = os.path.normpath("%s/%s/death_rates_rho=%g_p=%g_N=%g_K=%g.csv" % (curDir, folder, rho, p, N, K))
        r_v = np.loadtxt(path, delimiter = " ")
        avg_r += [np.average(r_v)]
    plt.plot(rho_vec[0:end], avg_r, "b")
    
    folder_2= "N_50_death_rates"
    end_2 = end
    avg_r = []
    for rho in rho_vec[0:end_2]:
        path = os.path.normpath("%s/%s/pr_death_rates_rho=%g_p=%g_N=%g_K=%g.csv" % (curDir, folder_2, rho, p, N, K))
        r_v = np.loadtxt(path, delimiter = " ")
        avg_r += [np.average(r_v)]
    
    plt.xlabel("vaccination rate $\\rho$")
    plt.ylabel("bare infection rate $r$")
    plt.plot(rho_vec[0:end_2], avg_r, "r")
    plt.show()
    i += 1
#
#plt.xlabel("$\\rho$")
#plt.ylabel("$r$")
##plt.title("Death rate $r$ vs. $\\rho$ for $p$ = 0.1 (green) and $p$ = 1 (red)")
#plt.show()

