from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt

#let's think of time steps as days

S0 = 1000
I0 = 1000
y0, t0 = [S0, I0, 0, S0 + I0], 0


a = 1/8 #removal rate I -> R, average disease length 8 days

b = 1 #incidence rate, # of people infected in a day

d = 0 #infectivity of treated T
nu = 1/8 #removal rate of treated T -> R
f_I = 0.1 #fraction of infected people who recover
f_T = 0.6 #fraction of treated people who recover

R0 = 200

g = (b*S0 - a*R0) / R0 #diagnosis rate I -> T, Depends on what we want R0 to be
print(g)
#if we want nonstandard incidence
l = 1 #number of people that an average person infects in unit time (standard incidence)
def beta(N):
    a = 0
    C = l*N**a
    return C/N


def f(t, y, alpha, beta, gamma, delta, nu, f_I, f_T):
    S, I, T, N = y
    S_prime = -beta*S*(I + delta*T)
    I_prime = beta*S*(I + delta*T) - (alpha + gamma)*I
    T_prime = gamma*I-nu*T
    N_prime = -(1 - f_I)*alpha*I - (1 - f_T)*nu*T
    return [S_prime, I_prime, T_prime, N_prime]
    
r = ode(f).set_integrator('dopri5')

r.set_initial_value(y0, t0).set_f_params(a, b, g, d, nu, f_I, f_T)

t1 = 5
dt = 0.1

t = np.empty((t1/dt) + 1)
S = np.empty((t1/dt) + 1)
I = np.empty((t1/dt) + 1)
T = np.empty((t1/dt) + 1)
N = np.empty((t1/dt) + 1)

i = 0
while r.successful() and r.t < t1:
    r.integrate(r.t+dt)
    #print("%g %s" % (r.t, r.y))
    t[i] = r.t
    S[i] = r.y[0]
    I[i] = r.y[1]
    T[i] = r.y[2]
    N[i] = r.y[3]
    i += 1

plt.close()
plt.plot(t, I, 'r')
plt.plot(t, S, 'b')
plt.plot(t, T, 'g')
plt.plot(t, N, 'k')
plt.xlabel("$t$")
plt.ylabel("$N$")
plt.title("Number vs Time of Alive (black), S (blue), I (red), T (green)")
#plt.savefig("S0_%s_I0_%s_a_%s_g_%s_d_%s_nu_%s_l_%s.png" % (S0, I0, a, g, d, nu, l))
#plt.savefig("./documents/Canopy/firstrun.png")
plt.show()
