from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt

def k_fate(ename, k_fate_0, k_fate_x, KD_Kx_fate, effector_cell_obs):
    return Expression(ename, (k_fate_0*KD_Kx_fate + k_fate_x*effector_cell_obs) / (KD_Kx_fate + effector_cell_obs))

##### NE #####

k_NE_div_0 = [1.0, 1.0] # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_NE_div_x = [2.0, 2.0]
KD_Kx_NE_div = [1000.0, 1000.0]
k_NE_die_0 = [0.9, 0.9]
k_NE_die_x = [0.1, 0.1]
KD_Kx_NE_die = [1000.0, 1000.0]

##### NEv1 #####

# Epithelium #
k_NEv1_div_0 = [1.0, 1.0] # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_NEv1_div_x = [2.0, 2.0]
KD_Kx_NEv1_div = [1000.0, 1000.0]
k_NEv1_die_0 = [0.9, 0.9]
k_NEv1_die_x = [0.1, 0.1]
KD_Kx_NEv1_die = [1000.0, 1000.0]

##### NEv2 #####

# Epithelium #
k_NEv2_div_0 = [1.0, 1.0] # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_NEv2_div_x = [2.0, 2.0]
KD_Kx_NEv2_div = [1000.0, 1000.0]
k_NEv2_die_0 = [0.9, 0.9]
k_NEv2_die_x = [0.1, 0.1]
KD_Kx_NEv2_die = [1000.0, 1000.0]

##### nonNE #####

# Epithelium #
k_nonNE_div_0 = [1.1, 1.1]  # [epithelium, stroma]
k_nonNE_div_x = [0.9, 0.9]
KD_Kx_nonNE_div = [1000.0, 1000.0]
k_nonNe_die = [0.1, 0.1]

##### NE <> NEv1 #####

# Epithelium #
kf_diff_ne_nev1 = [0.1, 0.1]  # [epithelium, stroma]
kr_diff_ne_nev1 = [0.1, 0.1]

##### NE <> NEv2 #####

# Epithelium #
kf_diff_ne_nev2 = [0.1, 0.1] # [epithelium, stroma]
kr_diff_ne_nev2 = [0.075, 0.075]

##### NEv1 <> NEv2 #####

# Epithelium #
kf_diff_nev1_nev2 = [0.1, 0.1] # [epithelium, stroma]
kr_diff_nev1_nev2 = [0.1, 0.1]

##### NEv1 >> nonNE #####

# Epithelium #
kf_diff_nev1_nonNe = [5.0, 5.0] # [epithelium, stroma]

Model()

mon = model.monomers
cmp = model.compartments
obs = model.observables
par = model.parameters
exp = model.expressions

Monomer('NE')
Monomer('NEv1')
Monomer('NEv2')
Monomer('NonNE')

Compartment('E', parent=None, dimension=3) # epithelium
Compartment('S', parent=None, dimension=3) # stroma

Parameter('NE_init_E', 100)
Initial(NE()**E, NE_init_E)

Parameter('NE_init_S', 0)
Initial(NE()**S, NE_init_S)

Observable('NE_obs_TOT', NE())
Observable('NEv1_obs_TOT', NEv1())
Observable('NEv2_obs_TOT', NEv2())
Observable('NonNE_obs_TOT', NonNE())
Observable('NE_all_TOT', NE()+NEv1()+NEv2())

[Observable('NE_obs_%s' % C.name, NE()**C) for C in cmp]
[Observable('NEv1_obs_%s' % C.name, NEv1()**C) for C in cmp]
[Observable('NEv2_obs_%s' % C.name, NEv2()**C) for C in cmp]
[Observable('NonNE_obs_%s' % C.name, NonNE()**C) for C in cmp]
[Observable('NE_all_%s' % C.name, NE()**C+NEv1()**C+NEv2()**C) for C in cmp]

# Parameter('k_ne_div', 1) 
[Parameter('k_NE_div_0_%s' % C.name, k_NE_div_0[i]) for i,C in enumerate(cmp)] 
[Parameter('k_NE_div_x_%s' % C.name, k_NE_div_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_NE_div_%s' % C.name, KD_Kx_NE_div[i]) for i,C in enumerate(cmp)]
[k_fate('k_NE_div_%s' % C.name, par['k_NE_div_0_%s' % C.name], par['k_NE_div_x_%s' % C.name], 
        par['KD_Kx_NE_div_%s' % C.name], obs['NonNE_obs_%s' % C.name]) for i,C in enumerate(cmp)] 
[Rule('NE_div_%s' % C.name, NE()**C >> NE()**C + NE()**C, exp['k_NE_div_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_NE_die_0_%s' % C.name, k_NE_die_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_NE_die_x_%s' % C.name, k_NE_die_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_NE_die_%s' % C.name, KD_Kx_NE_die[i]) for i,C in enumerate(cmp)]
[k_fate('k_NE_die_%s' % C.name, par['k_NE_die_0_%s' % C.name], par['k_NE_die_x_%s' % C.name], 
        par['KD_Kx_NE_die_%s' % C.name], obs['NonNE_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('NE_die_%s' % C.name, NE()**C >> None, exp['k_NE_die_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_NEv1_div_0_%s' % C.name, k_NEv1_div_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_NEv1_div_x_%s' % C.name, k_NEv1_div_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_NEv1_div_%s' % C.name, KD_Kx_NEv1_div[i]) for i,C in enumerate(cmp)]
[k_fate('k_NEv1_div_%s' % C.name, par['k_NEv1_div_0_%s' % C.name], par['k_NEv1_div_x_%s' % C.name], 
        par['KD_Kx_NEv1_div_%s' % C.name], obs['NonNE_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('NEv1_div_%s' % C.name, NEv1()**C >> NEv1()**C + NEv1()**C, exp['k_NEv1_div_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_NEv1_die_0_%s' % C.name, k_NEv1_die_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_NEv1_die_x_%s' % C.name, k_NEv1_die_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_NEv1_die_%s' % C.name, KD_Kx_NEv1_die[i]) for i,C in enumerate(cmp)]
[k_fate('k_NEv1_die_%s' % C.name, par['k_NEv1_die_0_%s' % C.name], par['k_NEv1_die_x_%s' % C.name], 
        par['KD_Kx_NEv1_die_%s' % C.name], obs['NonNE_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('NEv1_die_%s' % C.name, NEv1()**C >> None, exp['k_NEv1_die_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_NEv2_div_0_%s' % C.name, k_NEv2_div_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_NEv2_div_x_%s' % C.name, k_NEv2_div_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_NEv2_div_%s' % C.name, KD_Kx_NEv2_div[i]) for i,C in enumerate(cmp)]
[k_fate('k_NEv2_div_%s' % C.name, par['k_NEv2_div_0_%s' % C.name], par['k_NEv2_div_x_%s' % C.name], 
        par['KD_Kx_NEv2_div_%s' % C.name], obs['NonNE_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('NEv2_div_%s' % C.name, NEv2()**C >> NEv2()**C + NEv2()**C, exp['k_NEv2_div_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_NEv2_die_0_%s' % C.name, k_NEv2_die_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_NEv2_die_x_%s' % C.name, k_NEv2_die_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_NEv2_die_%s' % C.name, KD_Kx_NEv2_die[i]) for i,C in enumerate(cmp)]
[k_fate('k_NEv2_die_%s' % C.name, par['k_NEv2_die_0_%s' % C.name], par['k_NEv2_die_x_%s' % C.name], 
        par['KD_Kx_NEv2_die_%s' % C.name], obs['NonNE_obs_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('NEv2_die_%s' % C.name, NEv2()**C >> None, exp['k_NEv2_die_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_nonNE_div_0_%s' % C.name, k_nonNE_div_0[i]) for i,C in enumerate(cmp)]
[Parameter('k_nonNE_div_x_%s' % C.name, k_nonNE_div_x[i]) for i,C in enumerate(cmp)]
[Parameter('KD_Kx_nonNE_div_%s' % C.name, KD_Kx_nonNE_div[i]) for i,C in enumerate(cmp)]
[k_fate('k_nonNE_div_%s' % C.name, par['k_nonNE_div_0_%s' % C.name], par['k_nonNE_div_x_%s' % C.name], 
        par['KD_Kx_nonNE_div_%s' % C.name], obs['NE_all_%s' % C.name]) for i,C in enumerate(cmp)]
[Rule('NonNE_div_%s' % C.name, NonNE()**C >> NonNE()**C + NonNE()**C, exp['k_nonNE_div_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('k_nonNe_die_%s' % C.name, k_nonNe_die[i]) for i,C in enumerate(cmp)]
[Rule('NonNE_die_%s' % C.name, NonNE()**C >> None, par['k_nonNe_die_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('kf_diff_ne_nev1_%s' % C.name, kf_diff_ne_nev1[i]) for i,C in enumerate(cmp)]
[Parameter('kr_diff_ne_nev1_%s' % C.name, kr_diff_ne_nev1[i]) for i,C in enumerate(cmp)]
[Rule('NE_diff_NEv1_%s' % C.name, NE()**C | NEv1()**C, par['kf_diff_ne_nev1_%s' % C.name], 
      par['kr_diff_ne_nev1_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('kf_diff_ne_nev2_%s' % C.name, kf_diff_ne_nev2[i]) for i,C in enumerate(cmp)]
[Parameter('kr_diff_ne_nev2_%s' % C.name, kr_diff_ne_nev2[i]) for i,C in enumerate(cmp)]
[Rule('NE_diff_NEv2_%s' % C.name, NE()**C | NEv2()**C, par['kf_diff_ne_nev2_%s' % C.name], 
      par['kr_diff_ne_nev2_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('kf_diff_nev1_nev2_%s' % C.name, kf_diff_nev1_nev2[i]) for i,C in enumerate(cmp)]
[Parameter('kr_diff_nev1_nev2_%s' % C.name, kr_diff_nev1_nev2[i]) for i,C in enumerate(cmp)]
[Rule('NEv1_diff_NEv2_%s' % C.name, NEv1()**C | NEv2()**C, par['kf_diff_nev1_nev2_%s' % C.name], 
      par['kr_diff_nev1_nev2_%s' % C.name]) for i,C in enumerate(cmp)]

[Parameter('kf_diff_nev1_nonNe_%s' % C.name, kf_diff_nev1_nonNe[i]) for i,C in enumerate(cmp)]
[Rule('NEv1_diff_NonNE_%s' % C.name, NEv1()**C >> NonNE()**C, par['kf_diff_nev1_nonNe_%s' % C.name]) for i,C in enumerate(cmp)]

## Invasion (epithelium to stroma, one-way)

Parameter('k_NE_epi_to_stroma', 0.5)
Parameter('k_NEv1_epi_to_stroma', 0.5)
Parameter('k_NEv2_epi_to_stroma', 0.5)
Parameter('k_nonNE_epi_to_stroma', 0.5)

Rule('NE_epi_to_stroma', NE()**E >> NE()**S, k_NE_epi_to_stroma)
Rule('NEv1_epi_to_stroma', NEv1()**E >> NEv1()**S, k_NEv1_epi_to_stroma)
Rule('NEv2_epi_to_stroma', NEv2()**E >> NEv2()**S, k_NEv2_epi_to_stroma)
Rule('nonNE_epi_to_stroma', NonNE()**E >> NonNE()**S, k_nonNE_epi_to_stroma)

tspan = np.linspace(0, 20, 101)

sim = ScipyOdeSimulator(model, verbose=True)
x = sim.run(tspan)

for C_name in [c.name for c in cmp] + ['TOT']:
    
    obs_name = ['NE_obs_%s' % C_name, 'NEv1_obs_%s' % C_name, 'NEv2_obs_%s' % C_name, 'NonNE_obs_%s' % C_name]

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
    







