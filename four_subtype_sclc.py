from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt

def k_fate(ename, k_fate_0, k_fate_x, KD_Kx_fate, effector_cell_obs):
    return Expression(ename, (k_fate_0*KD_Kx_fate + k_fate_x*effector_cell_obs) / (KD_Kx_fate + effector_cell_obs))

Model()

Monomer('A') # NE (ASCL1)
Monomer('N') # NEv1 (NEUROD1)
Monomer('A2') # NEv2 (ASCL1)
Monomer('Y') # NonNE (YAP1)

Parameter('A_init', 100)
Initial(A(), A_init)

Observable('A_obs', A())
Observable('N_obs', N())
Observable('A2_obs', A2())
Observable('Y_obs', Y())
Observable('NE_all', A()+N()+A2())

# Parameter('k_ne_div', 1) #division role for the N subtype
Parameter('k_A_div_0', 1) # TPCs divide approximately once per day in culture
Parameter('k_A_div_x', 2)
Parameter('KD_Kx_A_div', 1000)
k_fate('k_A_div', k_A_div_0, k_A_div_x, KD_Kx_A_div, Y_obs)
Rule('A_div', A() >> A() + A(), k_A_div)   # >> division

# Parameter('k_ne_die', 0.9)
Parameter('k_A_die_0', 0.9)
Parameter('k_A_die_x', 0.1)
Parameter('KD_Kx_A_die', 1000)
k_fate('k_A_die', k_A_die_0, k_A_die_x, KD_Kx_A_die, Y_obs)
Rule('A_die', A() >> None, k_A_die)  # death

# Parameter('k_N_div', 1)
Parameter('k_N_div_0', 1) # TPCs divide approximately once per day in culture
Parameter('k_N_div_x', 2)
Parameter('KD_Kx_N_div', 1000)
k_fate('k_N_div', k_N_div_0, k_N_div_x, KD_Kx_N_div, Y_obs)
Rule('N_div', N() >> N() + N(), k_N_div)

# Parameter('k_nev1_die', 0.9)
Parameter('k_N_die_0', 0.9)
Parameter('k_N_die_x', 0.1)
Parameter('KD_Kx_N_die', 1000)
k_fate('k_N_die', k_N_die_0, k_N_die_x, KD_Kx_N_die, Y_obs)
Rule('N_die', N() >> None, k_N_die)

# Parameter('k_nev2_div', 1)
Parameter('k_A2_div_0', 1) # TPCs divide approximately once per day in culture
Parameter('k_A2_div_x', 2)
Parameter('KD_Kx_A2_div', 1000)
k_fate('k_A2_div', k_A2_div_0, k_A2_div_x, KD_Kx_A2_div, Y_obs)
Rule('A2_div', A2() >> A2() + A2(), k_A2_div)

# Parameter('k_nev2_die', 0.9)
Parameter('k_A2_die_0', 0.9)
Parameter('k_A2_die_x', 0.1)
Parameter('KD_Kx_A2_die', 1000)
k_fate('k_A2_die', k_A2_die_0, k_A2_die_x, KD_Kx_A2_die, Y_obs)
Rule('A2_die', A2() >> None, k_A2_die)

# Parameter('k_nonNe_div', 0.9)
Parameter('k_Y_div_0', 1.1)
Parameter('k_Y_div_x', 0.9)
Parameter('KD_Kx_Y_div', 1000)
k_fate('k_Y_div', k_Y_div_0, k_Y_div_x, KD_Kx_Y_div, NE_all)
Rule('Y_div', Y() >> Y() + Y(), k_Y_div)

Parameter('k_Y_die', 0.1)
Rule('Y_die', Y() >> None, k_Y_die)

Parameter('kf_diff_A_N', 0.1)  # differentiation
Parameter('kr_diff_A_N', 0.1)
Rule('A_diff_N', A() | N(), kf_diff_A_N, kr_diff_A_N)   # A can turn into N and vice-versa

Parameter('kf_diff_A_A2', 0.1)
Parameter('kr_diff_A_A2', 0.075)
Rule('A_diff_A2', A() | A2(), kf_diff_A_A2, kr_diff_A_A2)

Parameter('kf_diff_N_A2', 0.1)
Parameter('kr_diff_N_A2', 0.1)
Rule('N_diff_A2', N() | A2(), kf_diff_N_A2, kr_diff_N_A2)

Parameter('kf_diff_N_Y', 5)
Rule('N_diff_Y', N() >> Y(), kf_diff_N_Y)      # N turns to Y but Y does not turn to N

tspan = np.linspace(0, 20, 101)

sim = ScipyOdeSimulator(model, tspan, verbose=True)
x = sim.run()

plt.figure()
for obs in model.observables[:4]:
    label = obs.name[:obs.name.find('_')]
    plt.plot(tspan, x.all[obs.name], lw=3, label=label)
plt.xlabel('time (d)', fontsize=16)
plt.ylabel('cell count', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend(loc=0)
plt.tight_layout()

# cell_tot = np.array([sum(x.observables[i]) for i in range(len(x.observables))])
cell_tot = sum(x.all[obs.name] for obs in [A_obs, N_obs, A2_obs, Y_obs])

plt.figure()
label = [obs.name[:obs.name.find('_')] for obs in model.observables[:4]]
plt.fill_between(tspan, x.all[model.observables[0].name] / cell_tot, label=label[0])
sum_prev = x.all[model.observables[0].name]
for i in range(1,len(model.observables[:4])-1):
    plt.fill_between(tspan, (x.all[model.observables[i].name] + sum_prev) / cell_tot, sum_prev / cell_tot, label=label[i])
    sum_prev += x.all[model.observables[i].name]
plt.fill_between(tspan, [1]*len(tspan), sum_prev / cell_tot, label=label[-1])
plt.xlabel('time (d)', fontsize=16)
plt.ylabel('cell fraction', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc=(0.75,0.6), framealpha=1)
plt.tight_layout()

plt.show()
    







