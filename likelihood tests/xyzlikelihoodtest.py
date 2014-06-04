import numpy as np
import pylab
import matplotlib.pyplot as plt



N=1000000

def loglikelihood(thisx, thisy, thisz):
    thisr2 = thisx*thisx + thisy*thisy + thisz*thisz 
    if thisr2 < 1:
        return 1
        #return thisr2 * np.sin(np.arccos(thisz/np.sqrt(thisr2)))
    else:
        return -np.inf

def nextstep(oldx,oldy,oldz):
    newx = oldx + (-1 + 2 * np.random.rand())*.1
    newy = oldy + (-1 + 2 * np.random.rand())*.1
    newz = oldz + (-1 + 2 * np.random.rand())*.1
    return newx, newy, newz
    
x=[.5]
y=[0]
z=[0]
loglikelihoods=[0]
r2=[.25]
r=[.5]
r3=[.125]
costheta=[0]

for n in range(N):
    testx, testy, testz = nextstep(x[n], y[n], z[n])
    testloglikelihood = loglikelihood(testx, testy, testz)
    checkloglikelihood = np.log(np.random.rand())
    if testloglikelihood - loglikelihoods[n] > checkloglikelihood:
        x.append(testx)
        y.append(testy)
        z.append(testz)
        loglikelihoods.append(testloglikelihood)
        r2.append(testx*testx + testy*testy + testz*testz)
    else:
        x.append(x[n])
        y.append(y[n])
        z.append(z[n])
        loglikelihoods.append(loglikelihoods[n])
        r2.append(r2[n])
    r.append(np.sqrt(r2[n]))
    r3.append(r[-1]*r2[-1])
    costheta.append(z[-1]/r[-1])
pylab.figure(0)
xd=np.arange(-1,1,.01)
xr=0.75*(1-np.power(xd,2))
n, bins, patches = plt.hist([x[N/2:N],y[N/2:N],z[N/2:N]],np.arange(-1,1,.1),normed=1)
plt.title("X, Y, Z")
plt.plot(xd,xr)
pylab.figure(1)
plt.title("r^3")
n, bins, patches = plt.hist(r3[N/2:N],np.arange(0,1,.02),normed=1)
pylab.figure(2)
plt.title("cos(theta)")
n, bins, patches = plt.hist(costheta[N/2:N],np.arange(-1.05,1.05,.05),normed=1)

plt.show()
