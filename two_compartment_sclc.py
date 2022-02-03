from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt

def k_fate(ename, k_fate_0, k_fate_x, KD_Kx_fate, effector_cell_obs):
    return Expression(ename, (k_fate_0*KD_Kx_fate + k_fate_x*effector_cell_obs) / (KD_Kx_fate + effector_cell_obs))

##### NE ##### (A)

k_A_div_0 = [1.0, 1.0]  # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_A_div_x = [2.0, 2.0]
KD_Kx_A_div = [1000.0, 1000.0]
k_A_die_0 = [0.9, 0.9]
k_A_die_x = [0.1, 0.1]
KD_Kx_A_die = [1000.0, 1000.0]

##### NEv1 ##### (N)

k_N_div_0 = [1.0, 1.0]  # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_N_div_x = [2.0, 2.0]
KD_Kx_N_div = [1000.0, 1000.0]
k_N_die_0 = [0.9, 0.9]
k_N_die_x = [0.1, 0.1]
KD_Kx_N_die = [1000.0, 1000.0]

##### NEv2 ##### (A2)

k_A2_div_0 = [1.0, 1.0]  # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_A2_div_x = [2.0, 2.0]
KD_Kx_A2_div = [1000.0, 1000.0]
k_A2_die_0 = [0.9, 0.9]
k_A2_die_x = [0.1, 0.1]
KD_Kx_A2_die = [1000.0, 1000.0]

##### nonNE ##### (Y)

k_Y_div_0 = [1.1, 1.1]  # [epithelium, stroma]
k_Y_div_x = [0.9, 0.9]
KD_Kx_Y_div = [1000.0, 1000.0]
k_Y_die = [0.1, 0.1]

##### NE <> NEv1 #####

kf_diff_A_N = [0.1, 0.1]  # [epithelium, stroma]
kr_diff_A_N = [0.1, 0.1]

##### NE <> NEv2 #####

kf_diff_A_A2 = [0.1, 0.1] # [epithelium, stroma]
kr_diff_A_A2 = [0.075, 0.075]

##### NEv1 <> NEv2 #####

kf_diff_N_A2 = [0.1, 0.1] # [epithelium, stroma]
kr_diff_N_A2 = [0.1, 0.1]

##### NEv1 >> nonNE #####

kf_diff_N_Y = [5.0, 5.0] # [epithelium, stroma]

Model()

mon = model.monomers
cmp = model.compartments
obs = model.observables
par = model.parameters
exp = model.expressions

Monomer('A')
Monomer('N')
Monomer('A2')
Monomer('Y')

Compartment('E', parent=None, dimension=3)  # epithelium
Compartment('S', parent=None, dimension=3)  # stroma

Parameter('A_init_E', 100)
Initial(A()**E, A_init_E)

Parameter('A_init_S', 0)
Initial(A()**S, A_init_S)

Observable('A_obs_TOT', A())
Observable('N_obs_TOT', N())
Observable('A2_obs_TOT', A2())
Observable('Y_obs_TOT', Y())
Observable('NE_all_TOT', A()+N()+A2())

[Observable('A_obs_%s' % C.name, A()**C) for C in cmp]
[Observable('N_obs_%s' % C.name, N()**C) for C in cmp]
[Observable('A2_obs_%s' % C.name, A2()**C) for C in cmp]
[Observable('Y_obs_%s' % C.name, Y()**C) for C in cmp]
[Observable('NE_all_%s' % C.name, A()**C+N()**C+A2()**C) for C in cmp]


[Parameter('k_A_div_0_%s' % C.name, k_A_div_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_A_div_x_%s' % C.name, k_A_div_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_A_div_%s' % C.name, KD_Kx_A_div[i]) for i,C in enumerate(cmp)]
[k_fate('k_A_div_%s' % C.name, par['k_A_div_0_%s' % C.name], par['k_A_div_x_%s' % C.name],
        par['KD_Kx_A_div_%s' % C.name], obs['Y_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('A_div_%s' % C.name, A()**C >> A()**C + A()**C, exp['k_A_div_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_A_die_0_%s' % C.name, k_A_die_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_A_die_x_%s' % C.name, k_A_die_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_A_die_%s' % C.name, KD_Kx_A_die[i]) for i,C in enumerate(cmp)]
[k_fate('k_A_die_%s' % C.name, par['k_A_die_0_%s' % C.name], par['k_A_die_x_%s' % C.name],
        par['KD_Kx_A_die_%s' % C.name], obs['Y_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('A_die_%s' % C.name, A()**C >> None, exp['k_A_die_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_N_div_0_%s' % C.name, k_N_div_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_N_div_x_%s' % C.name, k_N_div_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_N_div_%s' % C.name, KD_Kx_N_div[i]) for i,C in enumerate(cmp)]
[k_fate('k_N_div_%s' % C.name, par['k_N_div_0_%s' % C.name], par['k_N_div_x_%s' % C.name],
        par['KD_Kx_N_div_%s' % C.name], obs['Y_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('N_div_%s' % C.name, N()**C >> N()**C + N()**C, exp['k_N_div_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_N_die_0_%s' % C.name, k_N_die_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_N_die_x_%s' % C.name, k_N_die_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_N_die_%s' % C.name, KD_Kx_N_die[i]) for i,C in enumerate(cmp)]
[k_fate('k_N_die_%s' % C.name, par['k_N_die_0_%s' % C.name], par['k_N_die_x_%s' % C.name],
        par['KD_Kx_N_die_%s' % C.name], obs['Y_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('N_die_%s' % C.name, N()**C >> None, exp['k_N_die_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_A2_div_0_%s' % C.name, k_A2_div_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_A2_div_x_%s' % C.name, k_A2_div_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_A2_div_%s' % C.name, KD_Kx_A2_div[i]) for i,C in enumerate(cmp)]
[k_fate('k_A2_div_%s' % C.name, par['k_A2_div_0_%s' % C.name], par['k_A2_div_x_%s' % C.name],
        par['KD_Kx_A2_div_%s' % C.name], obs['Y_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('A2_div_%s' % C.name, A2()**C >> A2()**C + A2()**C, exp['k_A2_div_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_A2_die_0_%s' % C.name, k_A2_die_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_A2_die_x_%s' % C.name, k_A2_die_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_A2_die_%s' % C.name, KD_Kx_A2_die[i]) for i,C in enumerate(cmp)]
[k_fate('k_A2_die_%s' % C.name, par['k_A2_die_0_%s' % C.name], par['k_A2_die_x_%s' % C.name],
        par['KD_Kx_A2_die_%s' % C.name], obs['Y_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('A2_die_%s' % C.name, A2()**C >> None, exp['k_A2_die_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_Y_div_0_%s' % C.name, k_Y_div_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_Y_div_x_%s' % C.name, k_Y_div_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_Y_div_%s' % C.name, KD_Kx_Y_div[i]) for i,C in enumerate(cmp)]
[k_fate('k_Y_div_%s' % C.name, par['k_Y_div_0_%s' % C.name], par['k_Y_div_x_%s' % C.name],
        par['KD_Kx_Y_div_%s' % C.name], obs['NE_all_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('Y_div_%s' % C.name, Y()**C >> Y()**C + Y()**C, exp['k_Y_div_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_Y_die_%s' % C.name, k_Y_die[i]) for i,C in enumerate(cmp)]
[Rule('Y_die_%s' % C.name, Y()**C >> None, par['k_Y_die_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('kf_diff_A_N_%s' % C.name, kf_diff_A_N[i]) for i,C in enumerate(cmp)]
[Parameter('kr_diff_A_N_%s' % C.name, kr_diff_A_N[i]) for i,C in enumerate(cmp)]
[Rule('A_diff_N_%s' % C.name, A()**C | N()**C, par['kf_diff_A_N_%s' % C.name],
      par['kr_diff_A_N_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('kf_diff_A_A2_%s' % C.name, kf_diff_A_A2[i]) for i,C in enumerate(cmp)]
[Parameter('kr_diff_A_A2_%s' % C.name, kr_diff_A_A2[i]) for i,C in enumerate(cmp)]
[Rule('A_diff_A2_%s' % C.name, A()**C | A2()**C, par['kf_diff_A_A2_%s' % C.name],
      par['kr_diff_A_A2_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('kf_diff_N_A2_%s' % C.name, kf_diff_N_A2[i]) for i,C in enumerate(cmp)]
[Parameter('kr_diff_N_A2_%s' % C.name, kr_diff_N_A2[i]) for i,C in enumerate(cmp)]
[Rule('N_diff_A2_%s' % C.name, N()**C | A2()**C, par['kf_diff_N_A2_%s' % C.name],
      par['kr_diff_N_A2_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('kf_diff_N_Y_%s' % C.name, kf_diff_N_Y[i]) for i,C in enumerate(cmp)]
[Rule('N_diff_Y_%s' % C.name, N()**C >> Y()**C, par['kf_diff_N_Y_%s' % C.name]) for i,C in enumerate(cmp)]

## Invasion (epithelium to stroma, one-way)

Parameter('k_A_epi_to_stroma', 0.5)
Parameter('k_N_epi_to_stroma', 0.5)
Parameter('k_A2_epi_to_stroma', 0.5)
Parameter('k_Y_epi_to_stroma', 0.5)

Rule('A_epi_to_stroma', A()**E >> A()**S, k_A_epi_to_stroma)
Rule('N_epi_to_stroma', N()**E >> N()**S, k_N_epi_to_stroma)
Rule('A2_epi_to_stroma', A2()**E >> A2()**S, k_A2_epi_to_stroma)
Rule('Y_epi_to_stroma', Y()**E >> Y()**S, k_Y_epi_to_stroma)

tspan = np.linspace(0, 20, 101)

sim = ScipyOdeSimulator(model, verbose=True)
x = sim.run(tspan)

for C_name in [c.name for c in cmp] + ['TOT']:
    
    obs_name = ['A_obs_%s' % C_name, 'N_obs_%s' % C_name, 'A2_obs_%s' % C_name, 'Y_obs_%s' % C_name]

    plt.figure()
    label = []
    for name in obs_name:
        label.append(name[:name.find('_')])
        if C_name == 'E':
            label[-1] += '_epith' 
        elif C_name == 'S':
            label[-1] += '_stroma'
        else:
            label[-1] += '_total'
        plt.plot(tspan, x.all[name], lw=3, label=label[-1])
    plt.xlabel('time (d)', fontsize=16)
    plt.ylabel('cell count', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc=0)
    plt.tight_layout()
    
    cell_tot = sum(x.all[name] for name in obs_name)
    
    plt.figure()
    plt.fill_between(tspan, x.all[obs_name[0]] / cell_tot, label=label[0])
    sum_prev = x.all[obs_name[0]]
    for i in range(1,len(obs_name)-1):
        plt.fill_between(tspan, (x.all[obs_name[i]] + sum_prev) / cell_tot, sum_prev / cell_tot, label=label[i])
        sum_prev += x.all[obs_name[i]]
    plt.fill_between(tspan, [1]*len(tspan), sum_prev / cell_tot, label=label[-1])
    plt.xlabel('time (d)', fontsize=16)
    plt.ylabel('cell fraction', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc=(0.75,0.6), framealpha=1)
    plt.tight_layout()

plt.show()
    



