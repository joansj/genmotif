import sys,os
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.cm as cm

if len(sys.argv)!=7:
	print 'python '+sys.argv[0]+' <stream_filename> <wmin> <wmax> <num_motifs> <freq_motifs> <max_time>'
	sys.exit()

fn_stream=sys.argv[1]
wmin=int(sys.argv[2])
wmax=int(sys.argv[3])
nclu=int(sys.argv[4])
freq=int(sys.argv[5])
tmax=sys.argv[6]
fn_out='output'
fn_bin='genmotif'
loglogfit=True
znorm=True
plot_average=True
maxplots=25

if '/' not in sys.argv[0]: path='./'
else:
	tmp=sys.argv[0].split('/')
	path=''
	for t in tmp[:-1]: path+=t+'/'
fn_out=path+fn_out
fn_bin=path+fn_bin
fn_conf=path+'parameters.conf'

print 'Cleaning...'
os.system('rm -f '+fn_out+'*.*')

print 'Load time series preview...'
stream=np.loadtxt(fn_stream,delimiter=',')
if len(stream.shape)==1: stream=np.array([stream]).T
print stream.shape
tslen,tsdim=stream.shape

print '\nRun...'
call=fn_bin+' -s '+fn_stream+' -l '+str(tslen)+' '+str(tsdim)+' -w '+str(wmin)+' '+str(wmax)+' -m '+str(nclu)+' '+str(freq)+' -t '+tmax+' -o '+fn_out+' -c '+fn_conf+' -v 1'
print call
os.system(call)

print '\nPlots...'
plt.figure()
plt.plot(stream[:np.min([10000,len(stream)])])

# Fitness curve
fitdata=np.loadtxt(fn_out+'.fit.txt',delimiter='\t')
fig,ax=plt.subplots(1)
fig.set_tight_layout(True)
plt.plot(fitdata[:,0],fitdata[:,1],drawstyle='steps-post')
if loglogfit:
	ax.set_xscale('log')
	ax.set_yscale('log')
plt.ylim(np.min(fitdata[:,1])*0.95,np.max(fitdata[:,1])*1.05)
plt.xlabel('Execution time [s]')
plt.ylabel('Objective function')

# Patterns

# Load data
result=np.loadtxt(fn_out+'.sol.txt',delimiter='\t')
X={}
for k in range(result.shape[0]):
	i=int(result[k,0])
	w=int(result[k,1])
	c=int(result[k,2])
	if c<0: continue
	if c>=maxplots: break
	if c not in X: X[c]=[]
	zz=[]
	for j in range(tsdim):
		finterp=interpolate.interp1d(np.linspace(0,1,w),stream[i:i+w,j])
		z=finterp(np.linspace(0,1,wmax))
		if znorm:
			z-=np.mean(z)
			sd=np.std(z)
			if sd>0: z/=float(sd)
		zz.append(z)
	X[c].append(np.array(zz).T)
uC=X.keys()
uC.sort()

# Plot
fig=plt.figure()
fig.set_tight_layout(True)
if znorm: plt.title('z-normalized patterns')
else: plt.title('unnormalized patterns')
nrows=int(np.ceil(np.sqrt(len(uC))))
ncols=int(np.ceil(float(len(uC))/float(nrows)))
colors=cm.jet(np.linspace(0,1,len(uC)*tsdim+1))[:-1,:]
for i in range(len(uC)):
	c=uC[i]
	plt.subplot(nrows,ncols,i+1)
	for t in X[c]:
		for d in range(t.shape[1]):
			plt.plot(t[:,d],color=colors[len(uC)*d+i,:],linewidth=1)
	if plot_average:
		avg=np.zeros((wmax,tsdim))
		for t in X[c]: avg+=t
		avg/=float(len(X[c]))
		for d in range(avg.shape[1]): plt.plot(avg[:,d],color='k',linewidth=2)
	if znorm: plt.ylim(-4,4)
	#plt.xlabel('('+str(len(X[c]))+' members)')
	plt.xlabel('Cluster '+str(c))

plt.show()

print ''
