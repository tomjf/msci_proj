import matplotlib.pyplot as plt

x=[1,2,3,4,5,6,7,8,9]
y=[10,20,30,40,50,60,70,80,90]

a = 1.0
b = 'sn'

plt.scatter(x,y)
plt.title("Jacobi Elliptic" + "__"+ b + "__" + "m = " + str(a))
plt.xlabel('Amplitude')
plt.ylabel('ns')
plt.axhline(y=0.96,color='r')
plt.savefig("Jacobi Elliptic" + "__"+ b + "__" + "m = " + str(a) + '.png')
