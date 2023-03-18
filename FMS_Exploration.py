import matplotlib.pyplot as plt
import numpy as np
import csv
from glob import glob
import peakutils
from scipy.optimize import curve_fit
from scipy.special import wofz,erf
from lmfit.models import Model
import time
import random
from scipy.special import jv

#Plotting stuff

# plt.figure(figsize=(12,9))
# plt.figure(figsize=(10,7.5))
plt.figure(figsize=(10,7.5))


ax=plt.subplot(111)
# ax.spines["top"].set_visible(False)
# ax.spines["bottom"].set_visible(False)
# ax.spines["right"].set_visible(False)
# ax.spines["left"].set_visible(False)

ax.get_xaxis().tick_bottom()    
ax.get_yaxis().tick_left()  


# This code is going to investigate, numerically, the 
# recovered beatnote spectrum for FMS in the limits
# that omega_mod is both greater and less than 
# the linewidth


def Lorentzian_Absorption(array,gamma,center):
	out=[np.e**(-0.5*gamma/((val-center)**2+0.25*gamma**2)) for val in array]
	out=normalizer(out)
	return out

def Lorentzian_Absorption_at_f(freq,gamma,center):
	out=np.e**(-0.5*gamma/((freq-center)**2+0.25*gamma**2))
	return out

def Lorentzian_Dispersion(array,gamma,center):
	out=[0.5*gamma*(val-center)/((val-center)**2+0.25*gamma**2) for val in array]
	out=normalizer(out)
	return out

def Lorentzian_Dispersion_at_f(freq,gamma,center):
	out=0.5*gamma*(freq-center)/((freq-center)**2+0.25*gamma**2)
	return out

# def Modulated_Laser(array,gamma,center,mod_freq,harmonic_order,M):
# 	light_field=np.zeros(len(array))
# 	for i in range(harmonic_order):
# 		if i>0:
# 			light_field=[val1+jv(i,M)*val2 for val1,val2 in zip(light_field,Lorentzian_Absorption(array,gamma,center+i*mod_freq))]
# 			light_field=[val1+jv(i,M)*((-1)**i)*val2 for val1,val2 in zip(light_field,Lorentzian_Absorption(array,gamma,center-i*mod_freq))]
# 		else:
# 			light_field=[val1+jv(i,M)*val2 for val1,val2 in zip(light_field,Lorentzian_Absorption(array,gamma,center))]
# 	return light_field

def Modulated_Laser(frequencies,gamma,center,mod_freq,harmonic_order,M):
	light_field=np.zeros(len(frequencies))
	if harmonic_order>=0:
		light_field[frequency_to_index(0+center)]=jv(0,M)
	if harmonic_order>=1:
		light_field[frequency_to_index(mod_freq+center)]=jv(1,M)
		light_field[frequency_to_index(-mod_freq+center)]=-jv(1,M)
	if harmonic_order>=2:
		light_field[frequency_to_index(2*mod_freq+center)]=jv(2,M)
		light_field[frequency_to_index(2*mod_freq+center)]=-jv(2,M)
	if harmonic_order>=3:
		light_field[frequency_to_index(3*mod_freq+center)]=jv(3,M)
		light_field[frequency_to_index(-3*mod_freq+center)]=-jv(3,M)
	if harmonic_order>=4:
		print('need to update modulated laser function')
	return light_field

def frequency_to_index(desired_frequency):
	return int((desired_frequency+f_span/2)*1/df)

def BeatNote_Demod(frequencies,array,filter,offset):
	# pos_freq_index=f.index(filter)
	# neg_freq_index=f.index(-filter)
	pos_freq_index=np.where((frequencies>filter-0.01+offset)&(frequencies<filter+.01+offset))[0][0]
	# print(pos_freq_index)
	neg_freq_index=np.where((frequencies>-filter-0.01+offset)&(frequencies<-filter+.01+offset))[0][0]
	# print(neg_freq_index)
	delta=3
	beat_spectrum=np.asarray(array)
	# print('rf power in positive lobe'+str(sum(beat_spectrum[pos_freq_index-delta:pos_freq_index+delta])))
	# print('rf power in negative lobe'+str(sum(beat_spectrum[neg_freq_index-delta:neg_freq_index+delta])))
	rf_pow=sum(beat_spectrum[pos_freq_index-delta:pos_freq_index+delta])+sum(beat_spectrum[neg_freq_index-delta:neg_freq_index+delta])
	return rf_pow

def normalizer(array):
	normed=[val/max(array) for val in array]
	return normed

def transmitted_pow(light,absorber):
	trans_pow=[val1*val2 for val1,val2 in zip(light,absorber)]
	# trans_pow=[val1*np.e**(-val2) for val1,val2 in zip(light,absorber)]
	return trans_pow

# def error_signal(frequencies,light_field_params,absorber,filter):
# 	error_sig=[]
# 	f3=np.linspace(-5,5,25)
# 	mod_index=1
# 	mod_freq=1
# 	demod_freq=mod_freq
# 	light_field=Modulated_Laser(f,gamma=0.01,center=0,mod_freq=mod_freq,harmonic_order=2,M=mod_index)
# 	for freq in f3:
# 		after_sample_shifted=transmitted_pow(Modulated_Laser(f,gamma=0.01,center=0+freq,mod_freq=mod_freq,harmonic_order=2,M=mod_index),absorption1)
# 		convolved_transmitted_pow=np.convolve(after_sample_shifted,after_sample_shifted,mode='same')
# 		error_sig.append(BeatNote_Demod(f,convolved_transmitted_pow,demod_freq,freq))
# 	plt.plot(f3,error_sig)
# 	plt.show()

def error_signal(frequencies,light_field_params,absorber,filter,mod_freq):
	error_sig=[]
	f3=np.linspace(-f_span/40,f_span/40,int(f_span/(df*20)))
	# f3=frequencies
	mod_index=0.1
	# mod_freq=1
	# demod_freq=mod_freq
	demod_freq=filter
	light_field=Modulated_Laser(frequencies,gamma=0.0,center=0,mod_freq=mod_freq,harmonic_order=1,M=mod_index)
	for freq in f3:
		# print(freq)
		absorption2=Lorentzian_Absorption(f,gamma=1,center=0+freq)
		after_sample=transmitted_pow(Modulated_Laser(f,gamma=0,center=0,mod_freq=mod_freq,harmonic_order=1,M=1),absorption2)
		convolved_laser=np.convolve(after_sample,after_sample,mode='same')
		error_sig.append(BeatNote_Demod(f,convolved_laser,filter=filter,offset=0))
	return f3,error_sig
	# plt.plot(f,after_sample)
	# plt.plot(f,convolved_laser)
	# plt.plot(f3,error_sig)
	# plt.plot(f,Lorentzian_Absorption(f,gamma=1,center=0))
	# plt.show()


df=0.1
f_span=100
x=np.linspace(0,f_span,int(f_span/df+1))
f=np.linspace(-f_span/2,f_span/2,int(f_span/df+1))

absorption1=Lorentzian_Absorption(f,gamma=1,center=0)
# absorption1=[np.e**(-val) for val in absorption1]
laser1=Modulated_Laser(f,gamma=0,center=0,mod_freq=1,harmonic_order=1,M=0.1)
# plt.plot(f,laser1)
# plt.show()

# after_sample=transmitted_pow(laser1,absorption1)
# plt.plot(f,after_sample)
# plt.plot(f,absorption1)
# plt.show()

# f_e1,error1=error_signal(f,0,absorption1,filter=1,mod_freq=1)
# f_e2,error2=error_signal(f,0,absorption1,filter=0.9,mod_freq=0.9)
# f_e3,error3=error_signal(f,0,absorption1,filter=0.8,mod_freq=0.8)
# f_e4,error1=error_signal(f,0,absorption1,filter=0.7,mod_freq=0.7)
# f_e5,error2=error_signal(f,0,absorption1,filter=0.6,mod_freq=0.6)
# f_e6,error3=error_signal(f,0,absorption1,filter=0.5,mod_freq=0.5)

# plt.plot(f,absorption1)


#check second harmonic!


# for i in range(13):
# 	f_0,error_0=error_signal(f,0,absorption1,filter=2*1.5-i/10,mod_freq=1.5-i/10)
# 	plt.plot(f_0,error_0,label='mod freq '+str(2*1.5-i/10))
# plt.xlim(-5,5)
# plt.legend()
# plt.show()

# convolved_laser1=np.convolve(after_sample,after_sample,mode='same')
# print(BeatNote_Demod(f,convolved_laser1,1,0))

# for i in range(5):
# 	after_sample=transmitted_pow(Modulated_Laser(f,gamma=0,center=0,mod_freq=1,harmonic_order=1,M=1),Lorentzian_Absorption(f,gamma=1,center=0+i/5))
# 	convolved_laser=np.convolve(after_sample,after_sample,mode='same')
# 	plt.plot(f,after_sample)
# 	plt.plot(f,convolved_laser)
# 	plt.show()
# 	print(i)
# 	print(BeatNote_Demod(f,convolved_laser,1,i))

#Analytic version from Hall 1985

def Hall_abs(frequencies,omega_mod,Gamma_0,detuning):
	recovered_sig=[]
	for freq in frequencies:
		value_from_first_resonance=Lorentzian_Absorption_at_f(freq+omega_mod,Gamma_0,0)-Lorentzian_Absorption_at_f(freq-omega_mod,Gamma_0,0)
		value_from_detuned_resonance=Lorentzian_Absorption_at_f(freq+omega_mod,Gamma_0,detuning)-Lorentzian_Absorption_at_f(freq-omega_mod,Gamma_0,detuning)
		recovered_sig.append(0.25*value_from_first_resonance+value_from_detuned_resonance)
	return recovered_sig

#proximity of other resonances?

def Hall_disp(frequencies,omega_mod,Gamma_0,detuning):
	recovered_sig=[]
	for freq in frequencies:
		value_from_first_resonance=Lorentzian_Dispersion_at_f(freq+omega_mod,Gamma_0,0)-2*Lorentzian_Dispersion_at_f(freq+0,Gamma_0,0)+Lorentzian_Dispersion_at_f(freq-omega_mod,Gamma_0,0)
		value_from_detuned_resonance=Lorentzian_Dispersion_at_f(freq+omega_mod,Gamma_0,detuning)-2*Lorentzian_Dispersion_at_f(freq+0,Gamma_0,detuning)+Lorentzian_Dispersion_at_f(freq-omega_mod,Gamma_0,detuning)
		recovered_sig.append(0.25*value_from_detuned_resonance+value_from_first_resonance)
	return recovered_sig

delta=-5

plt.plot(f,Hall_abs(f,1,1,delta),c='r',alpha=0.25)
plt.plot(f,Hall_abs(f,2,1,delta),c='r',alpha=0.25)
plt.plot(f,Hall_abs(f,3,1,delta),c='r',alpha=0.25)
plt.plot(f,Hall_abs(f,4,1,delta),c='r',alpha=0.25)
plt.plot(f,Hall_abs(f,5,1,delta),c='r',alpha=0.25)
# plt.plot(f,Hall_abs(f,6,1),c='r')
# plt.plot(f,Hall_abs(f,7,1),c='r')
# plt.plot(f,Hall_abs(f,8,1),c='r')
# plt.plot(f,Hall_abs(f,9,1),c='r')
plt.plot(f,Hall_abs(f,10,1,delta),c='r',alpha=0.25)
plt.plot(f,Hall_abs(f,15,1,delta),c='r',alpha=0.25)
plt.plot(f,Hall_abs(f,20,1,delta),c='r',alpha=0.25)

plt.plot(f,Hall_disp(f,1,1,delta))
plt.plot(f,Hall_disp(f,2,1,delta))
plt.plot(f,Hall_disp(f,3,1,delta))
plt.plot(f,Hall_disp(f,4,1,delta))
plt.plot(f,Hall_disp(f,5,1,delta))
# plt.plot(f,Hall_disp(f,6,1))
# plt.plot(f,Hall_disp(f,7,1))
# plt.plot(f,Hall_disp(f,8,1))
# plt.plot(f,Hall_disp(f,9,1))
plt.plot(f,Hall_disp(f,10,1,delta))
plt.plot(f,Hall_disp(f,15,1,delta))
plt.plot(f,Hall_disp(f,20,1,delta))
plt.show()

# plt.plot(f,absorption2,alpha=0.5,color='b')
# plt.plot(f,dispersion1,alpha=0.5,color='g')
# plt.plot(f,dispersion2,alpha=0.5,color='b')

# plt.plot(f,error_signal(f,0,absorption1,3))


# maximal_val=int(max(y3))
maximal_val=100










# plt.yticks(range(0, maximal_val, maximal_val//10), [str(y_val) for y_val in range(0, maximal_val, maximal_val//10)], fontsize=14)
plt.xticks(fontsize=14)
# for y in range(0, maximal_val, maximal_val//10):    
#     plt.plot(range(0, int(max(f))), [y] * len(range(0, int(max(f)))), "--", lw=0.5, color="black", alpha=0.3)  


# plt.plot(f,y3,alpha=0.5,color='r')

# plt.text(0.5*max(x)-0.5*min(x), max(y3), "Plot Title"    
#        "  ", fontsize=17, ha="center")

plt.ylabel("Amplitude", fontsize=14)
plt.xlabel("Frequency", fontsize=14)
plt.title("FMS Exploration", fontsize=18)
# ax.set_xlabel('X Label', fontsize=14)

# y1=np.asarray(y1)
# y2=np.asarray(y2)
# y3=np.asarray(y3)

# plt.fill_between(x, y3,y2, alpha=0.1,color="#f7847c")
# plt.fill_between(x, y2,y1, alpha=0.1,color="#7f83f0")


# plt.text(-100, 1, "Following Hall 1985")

plt.tight_layout(pad=1)

# plt.show()