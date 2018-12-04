import numpy as np
import pylab as plt

def tr(x,y):
    return y,-x
trans = np.vectorize(tr)

filename = "dense_matrix_mapping.txt"

with open( filename , 'r' ) as fo: 
    mg, ng = (int(v) for v in fo.readline().split()[1:] )
    mb, nb = (int(v) for v in fo.readline().split()[1:] )
    rp, cp = (int(v) for v in fo.readline().split()[1:] )

I, J, ip, jp, i, j = np.loadtxt( filename, int, unpack=True )
R = ip*cp + jp

plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

fig = plt.figure( figsize=(ng*0.6,mg*0.6))
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
ax.set_xlim( 0.5,ng+0.5 )
ax.set_ylim( -mg-0.5,-0.5 )

a = [v+1 for v in range(mg)]
b = [-v-1 for v in range(mg)]
ax.set_yticks( b )
ax.set_yticklabels( a  )

for a,b,c, il, jl in zip(I,J,R, i, j):
    d, e = trans(a,b)
    ax.text( d,e+0.1,c )
    ax.text( d,e-0.2, '({},{})'.format(il,jl) )

a, b = trans(I,J)
ax.scatter( a,b,c=R,cmap='hsv')
plt.show()
