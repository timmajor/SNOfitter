import numpy as np


N=100000
trial_x=np.random.rand(N)*2-1
trial_y=np.random.rand(N)*2-1
trial_z=np.random.rand(N)*2-1

sample_x=[]
sample_y=[]
sample_z=[]

for n in range(N):
    if (trial_x[n]*trial_x[n]+trial_y[n]*trial_y[n]+trial_z[n]*trial_z[n] < 1):
        sample_x.append(trial_x[n])
        sample_y.append(trial_y[n])
        sample_z.append(trial_z[n])

lnP_x=0
lnP_y=0
lnP_z=0
lnP_r=0
lnP_t=0
lnP_p=0
for i in range(len(sample_x)):
    lnP_x += np.log((1-sample_x[i]*sample_x[i])/np.pi*2)
    lnP_y += np.log((1-sample_y[i]*sample_y[i])/np.pi*2)
    lnP_z += np.log((1-sample_z[i]*sample_z[i])/np.pi*2)
    r = np.sqrt(sample_x[i]*sample_x[i]+sample_y[i]*sample_y[i]+sample_z[i]*sample_z[i])
    lnP_r += np.log(r*r)
    lnP_t += np.log(np.sin(np.arccos(sample_z[i]/r))/2)
    lnP_p += np.log(0.5/np.pi)

print "xyz", lnP_x+lnP_y+lnP_z
print "rtp", lnP_r+lnP_t+lnP_p
