import matplotlib.pyplot as plt
import numpy as np

a = [0.5,1.5,3.3,4.0]
x0 = 0.35
n = 100

#recursive function
def f(x,a):
    return a*x*(1-x)

#sensitivity to initial function
def sensitivity_map(a,x0,n):
    datax = np.zeros(n)
    datay = np.zeros_like(datax)
    for a in a:
        x = x0
        for i in range(n):
            xold = x
            x = f(x,a)
            datax[i] = i
            datay[i] = x
        plt.plot(datax,datay,label="x0=%f"%(a,))
    leg = plt.legend(loc='lower center', ncol=4, shadow=False, fancybox=True, bbox_to_anchor=(0.5,-0.3))
    leg.get_frame().set_alpha(0.5)
    plt.xlabel("$n$")
    plt.ylabel("$x_{n+1}$")
    plt.title("Sensitivity to Initial Condition")
    plt.savefig("sensitivity.eps", format="eps", bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.show()

#logistic map function
def logistic_map(a,x0,n):
    datax = np.zeros(n)
    datay = np.zeros_like(datax)
    for a in a:
        for i in range(n):
            x = i/n
            xold = x
            x = f(x,a)
            datax[i] = xold
            datay[i] = x
        plt.plot(datax,datay,label="a=%f"%(a,))
    leg = plt.legend(loc='lower center', ncol=4, shadow=False, fancybox=True, bbox_to_anchor=(0.5,-0.3))
    leg.get_frame().set_alpha(0.5)
    plt.xlabel("$x_{n}$")
    plt.ylabel("$x_{n+1}$")
    plt.title("Logistic Map")
    plt.savefig("logisticmap.eps", format="eps", bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.show()

#main
sensitivity_map(a,x0,n)
logistic_map(a,x0,n)