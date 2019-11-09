#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
#initializes necessary libraries

#%%
#This shell contains the basic function used throughout the program
def algebraic_relation(a,b,c,d,e,f):
    m = np.array([[1, 0, ((d*e-b*f)/(a*d-b*c))],[0, 1, ((a*f-c*e)/(a*d-b*c))]])
    pos_a = m[0,-1]
    pos_b = m[1,-1]
    return np.array([[pos_a],[pos_b]])
#algebraic_relation is used to find our a and b constants in data

def linear_fit(x, y, yerr):
    sum1 = np.sum(1/(yerr**2))
    sumx = np.sum(x/(yerr**2))
    sumy = np.sum(y/(yerr**2))
    sumxy = np.sum(x*y/(yerr**2))
    sumx2 = np.sum((x**2)/(yerr**2))
    delta = (sum1*sumx2) - ((sumx)**2)
    fit_erra = np.sqrt((1/delta)*sumx2)
    fit_errb = np.sqrt((1/delta)*sum1)
    dof = len(x) - 1
    
    return (algebraic_relation(sum1, sumx, sumx, sumx2, sumy, sumxy), np.array[fit_erra,fit_errb], dof)
#linear_fit finds the sums found from chi2, then prints out their algebraic relation,
#fitted error parameters, and degrees of freedom
def fitted_data(x, y, yerr):
    y_int = linear_fit(x, y, yerr)[0] #This is the value of a
    slope = linear_fit(x, y, yerr)[1] #This is the value of b
    
    plt.figure(1)
    plt.errorbar(x, y, yerr=yerr, xerr=None, fmt='bo', ecolor='r', label='Data w/ Error Bars')
    plt.plot(x, slope*x + y_int, label='Linear Fit')
    plt.title('Linear Fit on data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend
    plt.show()
    
    plt.figure(2)
    residual = ((y - y_int - (slope*x)) / (yerr))
    plt.plot(x, residual, 'o')
    plt.title('Residuals vs. x')
    plt.xlabel('x')
    plt.ylabel('Residuals')
    plt.show()
    
    chi2 = np.sum(((y - y_int - (slope*x))/(yerr))**2)
    dof = len(x)-1
    p = 1-stats.chi2.cdf(5.0**2,1)
    print('chi2 equals '+str(chi2))
    print('dof equals '+str(dof))
    print('probability level equals '+str(p))
    return 0
#fitted_data fits the data linearly, produces two plots (one of fitted data, one of residuals)
#and also prints the chi2, dof, and probability level of the data

#%%
#This shell contains linear fits for the TEC data
(V, A, err) = np.loadtxt('TEC_data.csv', delimiter = ',', unpack=True)

fitted_data(V, A, err)
y_int = linear_fit(V, A, err)[0] #This is the value of a
slope = linear_fit(V, A, err)[1] #This is the value of b

plt.figure(1)
plt.errorbar(V, A, yerr=err, xerr=None, fmt='bo', ecolor='r', label='Data w/ Error Bars')
plt.plot(V, slope*V + y_int, label='Linear Fit')
plt.title('Linear Fit on TEC data')
plt.xlabel('Voltage')
plt.ylabel('Current')
plt.legend
plt.show()

plt.figure(2)
residual = ((A - y_int - (slope*V)) / (err))
plt.plot(V, residual, 'o')
plt.title('Residuals vs. Voltage')
plt.xlabel('Voltage')
plt.ylabel('Residuals')
plt.show()

chi2 = np.sum(((A - y_int - (slope*V))/(err))**2)
dof = len(V) - 1
p = 1-stats.chi2.cdf(5.0**2,1)
print('chi2 equals '+str(chi2))
print('dof equals '+str(dof))
print('the condifence level is '+str(p))

#%%
#This shell contains linear fits for the Fl18 data
(t, N) = np.loadtxt('fluorine18.csv', delimiter = ',', unpack=True)
(err) = np.sqrt(N)

#from fluorine18 data, N0 is 1021 at time 0
Tau = (t/(np.log(1021/N)))

#linear equation for exponential decay is
#np.log(N)=np.log(1021)-(t/Tau)
#where y = np.log(N), x = t, y_int = np.log(1021), slope = -1/Tau
#y_int = np.log(1021)
#slope = -1/Tau

fitted_data(t, N, 1/err)
y_int = linear_fit(t, N, 1/err)[0] #This is the value of a
slope = linear_fit(t, N, 1/err)[1] #This is the value of b
#I used this instead of the commented y_int and slop above because I kept receiving a
#nan error when t = 0. (Dividing by 0 brings up the nan error)

plt.figure(1)
plt.errorbar(t, N, yerr=err, xerr=None, fmt='bo', ecolor='r', label='Data w/ Error Bars')
plt.plot(t, slope*t + y_int, label='Linear Fit')
plt.title('Linear Fit on Fl18 radioactive decay data')
plt.xlabel('Time')
plt.ylabel('Counts (N)')
plt.legend
plt.show()

plt.figure(2)
residual = ((N - y_int - (slope*t)) / (err))
plt.plot(t, residual, 'o')
plt.title('Residuals vs. Time')
plt.xlabel('Time')
plt.ylabel('Residuals')
plt.show()

chi2 = np.sum(((N - y_int - (slope*t))/(err))**2)
dof = len(N)-1
p = 1-stats.chi2.cdf(5.0**2,1)
print('chi2 equals '+str(chi2))
print('dof equals '+str(dof))
print('the confidence level is '+str(p))