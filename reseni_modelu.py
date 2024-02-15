import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
pocatecni_podminka = [1312, 1000, 13827, 2]	# [lisky, skunci, prasata, orli]
meze = [0,100]	# cas pocatku a konce simulace
n = 100 # pocet delicich bodu
def model(t, N,
	rf=0.32, 
	rs=0.32, 
	rp=0.78, 
	Kf=1544, 
	Ks=2490, 
	Kp=15189, 
	Beta_fs=0.36, 
	Beta_sf=2.76, 
	mu_f=0.086, 
	mu_s=0.159, 
	mu_p=0.019, 
	Phi=8.1, 
	sigma=3.1, 
	ny=0.09, 
	lambda_f=0.00077, 
	lambda_s=0.00025, 
	lambda_p=0.00077):
	[foxes, skunks, pigs, eagles]=N
	d_foxes = rf * foxes * (1- (foxes + Beta_fs * skunks) / Kf) - mu_f * (Phi * foxes / (Phi * foxes + sigma * skunks + pigs)) * eagles  * foxes
	d_skunks = rs * skunks * (1 - (skunks + Beta_sf * foxes) / Ks) - mu_s * (sigma * skunks / (Phi * foxes + sigma * skunks + pigs) ) * eagles * skunks
	d_pigs = rp * pigs * (1 - pigs / Kp) - mu_p *pigs / (Phi * foxes + sigma * skunks + pigs) * eagles * pigs
	d_eagles = eagles * ( lambda_f * mu_f * Phi * foxes * foxes + lambda_s * mu_s * sigma * skunks *skunks + lambda_p *mu_p *pigs * pigs ) / (Phi * foxes + sigma * skunks + pigs) - ny * eagles
	return[d_foxes, d_skunks, d_pigs, d_eagles]
t=np.linspace(*meze, n)
pocty = solve_ivp(model, meze, pocatecni_podminka, t_eval=t)
fig, ax = plt.subplots()
pocty.y[2]=pocty.y[2]/10
pocty.y[3]=pocty.y[3]*10
ax.plot(t,pocty.y.T, label=["lisky","skunci","prasata","orli"])
ax.legend()
ax.set_xlabel("cas")
ax.set_ylabel("pocet lisek, skunku, prasat/10, orlu*10")
plt.show()
