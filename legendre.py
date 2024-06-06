
yvals = []
yvals2 = []
xth = np.linspace(-1,1,41)

Pl_1 = legendre(1)
for j in range(Nth+1):
    yvals.append(special.lpmv(1,4,xth[j]))
 
 

print special.lpmv(1,2,1)
print special.lpmv(1,2,0.5)
print special.lpmv(1,2,0.8)
Pl_1 = legendre(2)
print Pl_1(1)
print Pl_1(0.5)
print Pl_1(0.8)
#plt.plot(xth,yvals)
#plt.show()