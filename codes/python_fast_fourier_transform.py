#!/usr/bin/python
# -*- coding: utf-8 -*- 

def packet(E0,t,w0,T):
    return E0*cos(w0*t)*exp(-t**2/T**2) 

def packet_spectrum(E0,w,w0,T):
    return (E0*T/(2.0*sqrt(2.0)))**2*(exp(-(w - w0)**2*T**2/4.0) + exp(-(w + w0)**2*T**2/4.0))**2

###############################################################################
# First we plot the analytical signal and spectral density
from pylab import *
rcParams['axes.grid'] = True
rcParams['font.size'] = 14
rcParams['axes.labelsize'] = 22

E0 = 1.0
w0 = 1.0
T = 10.0

###############################################################################
# time axis
n = 800
t = linspace(-150,150,n)
dt = t[1] - t[0]

###############################################################################
# angular frequency axis
dw = 2.0*pi/(t[-1] - t[0])
w = arange(-n/2,n/2)*dw

###############################################################################
signal = packet(E0,t,w0,T) 
ana_spec = packet_spectrum(E0,w,w0,T)
trans = fftshift(fft(signal))*dt/sqrt(2.0*pi)
num_spec = abs(trans)**2
inv_trans = ifft(ifftshift(trans))*dw/sqrt(2.0*pi)*n

###############################################################################
figure(figsize=figaspect(0.5))

# time series
subplot(121)
plot(t, signal)
plot(t, inv_trans, "ro")
xlabel(r"$t$") 
ylabel(r"$E(t)$")
xlim(-30,30)

# spectrum
subplot(122)
plot(w, ana_spec,'-',label='analytical')
plot(w, num_spec,'ro',label='numerical')
xlim(0,2)
xlabel(r"$\omega$") 
ylabel(r"$|E(\omega)|^2$")
legend(loc='best')

tight_layout()
show()
###############################################################################