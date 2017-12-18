#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 15:11:12 2017

@author: xyh
"""

from pylab import *
from scipy.integrate import romberg

clight = 2.9979e5;
H0=70;
Om=0.3;
Or=2.46e-5;
Ol=1.-Om-Or;

def aH_inv(lna):
    a = exp(lna)
    return 1./H0/(Om/a+Ol*a**2+Or/a**2)**0.5

def a2H_inv(a):
    return 1./H0/(Om*a+Ol*a**4+Or)**0.5

zmin=0.001
zmax=5.0
#lnz = linspace(log(zmin),log(zmax),100)
#z = exp(lnz)
z = linspace(zmin,zmax,20)

a = 1/(1+z)
lna = log(a)

dc1 = []
dc2 = []
mu1 = []
mu2 = []

for i in range(len(a)):
#    dc1_temp = romberg(a2H_inv,a[i],1.0,tol=1e-10,rtol=1e-10,divmax=100)
#    dc2_temp = romberg(aH_inv,lna[i],0.0,tol=1e-10,rtol=1e-10,divmax=100)
    dc1_temp = romberg(a2H_inv,a[i],1.0)
    dc2_temp = romberg(aH_inv,lna[i],0.0)
    dc1.append(dc1_temp)
    dc2.append(dc2_temp)

    mu1.append(5.*log10((1+z[i])*dc1_temp*clight)+25)
    mu2.append(5.*log10((1+z[i])*dc2_temp*clight)+25)
    
dc1 = array(dc1)
dc2 = array(dc2)
mu1 = array(mu1)
mu2 = array(mu2)

ddc = (dc1-dc2)
dmu = mu1-mu2

#semilogy(z,abs(ddc),label=r'abs-err')
#plot(z,ddc/dc1,label=r'rel-err1')
#semilogy(z,abs(ddc)/dc2,label=r'rel-err2')

figure(figsize=(14,6))

subplot(1,2,1)
plot(z,dc1,'x',label=r'dc1')
plot(z,dc2,'+',label=r'dc2')

subplot(1,2,2)
#semilogy(z,abs(dmu),label=r'mu-diff')
#semilogx(z,dmu,'-')
#loglog(z,abs(dmu),'-')
#plot(z,mu1,'.',label=r'mu1')
#plot(z,mu2,'.',label=r'mu2')
semilogx(z,mu1,'x',label=r'mu1')
semilogx(z,mu2,'+',label=r'mu2')

legend()
show()

#z = 1500
#a = 1./(1.+z)
#dc1_temp = romberg(a2H_inv,a,1.0,divmax=100)
#dc2_temp = romberg(aH_inv,log(a),0.0,divmax=100)
#
#print 'dc1 = %g,\ndc2 = %g\n'%(dc1_temp,dc2_temp)
#print 'ddc_ref = %g'%((dc1_temp-dc2_temp)/dc1_temp)