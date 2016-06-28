import sys,os
import numpy as np

if len(sys.argv)!=8:
    print 'python '+sys.argv[0]+' <input_time_series> <wmin> <wmax> <k> <f> <time> <output_basename>'
    sys.exit()

fn_in=sys.argv[1]
wmin=int(sys.argv[2])
wmax=int(sys.argv[3])
k=int(sys.argv[4])
f=int(sys.argv[5])
time=sys.argv[6]
fn_base_out=sys.argv[7]
fn_bin='./genmotif'
fn_conf='parameters.conf'
numruns=5

print 'Load preview...'
x=np.loadtxt(fn_in,delimiter=',')
if len(x.shape)==1: x=np.array([x]).T
print x.shape

print 'Run...'
for r in range(numruns):
    call=fn_bin+' -s '+fn_in+' -l '+str(x.shape[0])+' '+str(x.shape[1])+' -w '+str(wmin)+' '+str(wmax)+' -m '+str(k)+' '+str(f)+' -t '+time+' -o '+fn_base_out+'_'+str(r)+' -c '+fn_conf
    print r+1,'-->',call
    os.system(call)
print 'Done'
