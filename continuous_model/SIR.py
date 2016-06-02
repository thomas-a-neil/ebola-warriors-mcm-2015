from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt

S0 = 1000
I0 = 100

y0, t0 = [S0, I0], 0

def f(t, y, beta, alpha):
    S, I = y
    S_prime = -beta*S*I
    I_prime = beta*S*I - alpha*I
    return [S_prime, I_prime]
    
r = ode(f).set_integrator('dopri5')

R = 5
alpha = 1.0/8
beta = (R * alpha) / S0


r.set_initial_value(y0, t0).set_f_params(beta, alpha)

t1 = 30
dt = 0.2

t = np.empty((t1/dt) + 1)
S = np.empty((t1/dt) + 1)
I = np.empty((t1/dt) + 1)

i = 0
while r.successful() and r.t < t1:
    r.integrate(r.t+dt)
    #print("%g %s" % (r.t, r.y))
    t[i] = r.t
    S[i] = r.y[0]
    I[i] = r.y[1]
    i += 1

plt.plot(t, I, 'r')
plt.plot(t, S, 'b')
plt.xlabel("$t$")
plt.ylabel("$N$")
#plt.savefig("SIR_R_%g_I0_%g.png" % (R, I0))

plt.show()