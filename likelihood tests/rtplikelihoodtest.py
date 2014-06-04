import numpy as np
import pylab
import matplotlib.pyplot as plt



N=100000

def loglikelihood(thisr, thist, thisp): 
    if thisr < 1 and thisr > 0 and thist < np.pi and thist > 0 and thisp < 2*np.pi and thisp> 0:
        return np.log(thisr*thisr) +np.log(np.sin(thist))
    else:
        return -np.inf

def nextstep(oldr,oldt,oldp):
    newr = oldr + (-1 + 2 * np.random.rand())*.2
    newt = oldt + (-1 + 2 * np.random.rand())*.6
    newp = oldp + (-1 + 2 * np.random.rand())*1.2
    return newr, newt, newp
    
loglikelihoods=[0]
r=[.5]
t=[.5*np.pi]
p=[np.pi]
r3=[.125]
costheta=[0]
x=[-0.5]
y=[0]
z=[0]


for n in range(N):
    testr, testt, testp = nextstep(r[n], t[n], p[n])
    testloglikelihood = loglikelihood(testr, testt, testp)
    checkloglikelihood = np.log(np.random.rand())
    if testloglikelihood - loglikelihoods[n] > checkloglikelihood:
        r.append(testr)
        t.append(testt)
        p.append(testp)
        loglikelihoods.append(testloglikelihood)
    else:
        r.append(r[n])
        t.append(t[n])
        p.append(p[n])
        loglikelihoods.append(loglikelihoods[n])
    r3.append(np.power(r[-1],3))
    costheta.append(np.cos(t[-1]))
    z.append(costheta[-1]*r[-1])
    x.append(np.cos(p[-1])*np.sin(t[-1])*r[-1])
    y.append(np.sin(p[-1])*np.sin(t[-1])*r[-1])
    
pylab.figure(0)
xd=np.arange(-1,1,.01)
xr=0.75*(1-np.power(xd,2))
n, bins, patches = plt.hist([x[N/2:N],y[N/2:N],z[N/2:N]],np.arange(-1,1,.05),normed=1)
plt.title("X, Y, Z")
plt.plot(xd,xr)
pylab.figure(1)
plt.title("r^3")
n, bins, patches = plt.hist(r3[N/2:N],np.arange(0,1,.01),normed=1)
pylab.figure(2)
plt.title("cos(theta)")
n, bins, patches = plt.hist(costheta[N/2:N],np.arange(-1,1,.03),normed=1)
pylab.figure(3)
n, bins, patches = plt.hist(p[N/2:N],np.arange(0,np.pi*2,.06),normed=1)

plt.show()
