from pysb import *
from pysb.simulator import ScipyOdeSimulator
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt
from sympy import Piecewise

# "For harvesting invasive lung tumor tissue, RPM mice were sacrified at 4-9 weeks post-adenoviral infection"
# Ireland et al. (2020) Cancer Cell, https://doi.org/10.1016/j.ccell.2020.05.001


def k_fate(ename, k_fate_0, k_fate_x, kd_kx_fate, effector_cell_obs):
    return Expression(ename, (k_fate_0*kd_kx_fate + k_fate_x*effector_cell_obs) / (kd_kx_fate + effector_cell_obs))


# #### A #####

k_A_div_0 = [2.0, 1.0]  # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_A_div_x = [8.0, 5.0]  # [50.0, 1.0]
KD_Kx_A_div = [1000.0, 1000.0]
k_A_die_0 = [0.2, 0.2]  # [2.0, 2.0]
k_A_die_x = [2.0, 2.0]
KD_Kx_A_die = [1000.0, 1000.0]

# #### N #####

k_N_div_0 = [2.0, 1.0]  # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_N_div_x = [9.0, 5.0]  # [50.0, 1.0]
KD_Kx_N_div = [1000.0, 1000.0]
k_N_die_0 = [1.9, 1.9]
k_N_die_x = [2.25, 2.25]
KD_Kx_N_die = [1000.0, 1000.0]

# #### A2 #####

k_A2_div_0 = [2.0, 1.0]  # [epithelium, stroma] # TPCs divide approximately once per day in culture
k_A2_div_x = [8.0, 5.0]  # [50.0, 1.0]
KD_Kx_A2_div = [1000.0, 1000.0]
k_A2_die_0 = [0.2, 0.2]
k_A2_die_x = [2.0, 2.0]
KD_Kx_A2_die = [1000.0, 1000.0]

# #### Y #####

k_Y_div_0 = [0.4, 4.0]  # [epithelium, stroma]
k_Y_div_x = [0.3, 3.0]  # [2.0, 1.0]
KD_Kx_Y_div = [1000.0, 1000.0]  # [500.0, 10000.0]
k_Y_die = [0.0, 2.0]

# #### A <> N #####

kf_diff_A_N = [0.1, 0.5]  # [0.01, 0.1]  # [epithelium, stroma]
kr_diff_A_N = [0.5, 0.01]  # [0.1, 0.1]

# #### A <> A2 #####

kf_diff_A_A2 = [0, 0]  # [0.05, 0.05]  # [epithelium, stroma]
kr_diff_A_A2 = [0, 0]  # [0.075, 0.075]

# #### N <> A2 #####

kf_diff_N_A2 = [0.5, 0.01]  # [0.001, 0.1]  # [epithelium, stroma]
kr_diff_N_A2 = [0.1, 0.5]  # [0.1, 0.1]

# #### N <> Y #####

kf_diff_N_Y = [0, 0.5]  # [0.001, 0.1]  # [epithelium, stroma]
kr_diff_N_Y = [0.5, 0]  # [epithelium, stroma]

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

Parameter('A_init_E', 1)  # 258)
Parameter('N_init_E', 0)  # 41)
Parameter('A2_init_E', 0)  # 89)
Parameter('Y_init_E', 0)  # 13)

Initial(A()**E, A_init_E)
Initial(N()**E, N_init_E)
Initial(A2()**E, A2_init_E)
Initial(Y()**E, Y_init_E)

Parameter('A_init_S', 0)
Parameter('N_init_S', 0)
Parameter('A2_init_S', 0)
Parameter('Y_init_S', 0)

Initial(A()**S, A_init_S)
Initial(N()**S, N_init_S)
Initial(A2()**S, A2_init_S)
Initial(Y()**S, Y_init_S)

Observable('A_obs_TOT', A())
Observable('A2_obs_TOT', A2())
Observable('N_obs_TOT', N())
Observable('Y_obs_TOT', Y())
Observable('NE_all_TOT', A()+N()+A2())
Observable('Cells_TOT', A()+N()+A2()+Y())

[Observable('A_obs_%s' % C.name, A()**C) for C in cmp]
[Observable('A2_obs_%s' % C.name, A2()**C) for C in cmp]
[Observable('N_obs_%s' % C.name, N()**C) for C in cmp]
[Observable('Y_obs_%s' % C.name, Y()**C) for C in cmp]
[Observable('NE_all_%s' % C.name, A()**C+N()**C+A2()**C) for C in cmp]
[Observable('Cells_%s' % C.name, A()**C+N()**C+A2()**C+Y()**C) for C in cmp]

# #### A div and death #####

[Parameter('k_A_div_0_%s' % C.name, k_A_div_0[i]) for i, C in enumerate(cmp)]
[Parameter('k_A_div_x_%s' % C.name, k_A_div_x[i]) for i, C in enumerate(cmp)]
[Parameter('KD_Kx_A_div_%s' % C.name, KD_Kx_A_div[i]) for i, C in enumerate(cmp)]
[k_fate('k_A_div_%s' % C.name, par['k_A_div_0_%s' % C.name], par['k_A_div_x_%s' % C.name],
        par['KD_Kx_A_div_%s' % C.name], obs['Y_obs_%s' % C.name]) for i, C in enumerate(cmp)]
[Rule('A_div_%s' % C.name, A()**C >> A()**C + A()**C, exp['k_A_div_%s' % C.name]) for i, C in enumerate(cmp)]

[Parameter('k_A_die_0_%s' % C.name, k_A_die_0[i]) for i, C in enumerate(cmp)]
[Parameter('k_A_die_x_%s' % C.name, k_A_die_x[i]) for i, C in enumerate(cmp)]
[Parameter('KD_Kx_A_die_%s' % C.name, KD_Kx_A_die[i]) for i, C in enumerate(cmp)]
[k_fate('k_A_die_%s' % C.name, par['k_A_die_0_%s' % C.name], par['k_A_die_x_%s' % C.name],
        par['KD_Kx_A_die_%s' % C.name], obs['Y_obs_%s' % C.name]) for i, C in enumerate(cmp)]
[Rule('A_die_%s' % C.name, A()**C >> None, exp['k_A_die_%s' % C.name]) for i, C in enumerate(cmp)]

# #### N div and death #####

[Parameter('k_N_div_0_%s' % C.name, k_N_div_0[i]) for i, C in enumerate(cmp)]
[Parameter('k_N_div_x_%s' % C.name, k_N_div_x[i]) for i, C in enumerate(cmp)]
[Parameter('KD_Kx_N_div_%s' % C.name, KD_Kx_N_div[i]) for i, C in enumerate(cmp)]
[k_fate('k_N_div_%s' % C.name, par['k_N_div_0_%s' % C.name], par['k_N_div_x_%s' % C.name],
        par['KD_Kx_N_div_%s' % C.name], obs['Y_obs_%s' % C.name]) for i, C in enumerate(cmp)]
[Rule('N_div_%s' % C.name, N()**C >> N()**C + N()**C, exp['k_N_div_%s' % C.name]) for i, C in enumerate(cmp)]

[Parameter('k_N_die_0_%s' % C.name, k_N_die_0[i]) for i, C in enumerate(cmp)]
[Parameter('k_N_die_x_%s' % C.name, k_N_die_x[i]) for i, C in enumerate(cmp)]
[Parameter('KD_Kx_N_die_%s' % C.name, KD_Kx_N_die[i]) for i, C in enumerate(cmp)]
[k_fate('k_N_die_%s' % C.name, par['k_N_die_0_%s' % C.name], par['k_N_die_x_%s' % C.name],
        par['KD_Kx_N_die_%s' % C.name], obs['Y_obs_%s' % C.name]) for i, C in enumerate(cmp)]
[Rule('N_die_%s' % C.name, N()**C >> None, exp['k_N_die_%s' % C.name]) for i, C in enumerate(cmp)]

# #### A2 div and death #####

[Parameter('k_A2_div_0_%s' % C.name, k_A2_div_0[i]) for i, C in enumerate(cmp)]
[Parameter('k_A2_div_x_%s' % C.name, k_A2_div_x[i]) for i, C in enumerate(cmp)]
[Parameter('KD_Kx_A2_div_%s' % C.name, KD_Kx_A2_div[i]) for i, C in enumerate(cmp)]
[k_fate('k_A2_div_%s' % C.name, par['k_A2_div_0_%s' % C.name], par['k_A2_div_x_%s' % C.name],
        par['KD_Kx_A2_div_%s' % C.name], obs['Y_obs_%s' % C.name]) for i, C in enumerate(cmp)]
[Rule('A2_div_%s' % C.name, A2()**C >> A2()**C + A2()**C, exp['k_A2_div_%s' % C.name]) for i, C in enumerate(cmp)]

[Parameter('k_A2_die_0_%s' % C.name, k_A2_die_0[i]) for i, C in enumerate(cmp)]
[Parameter('k_A2_die_x_%s' % C.name, k_A2_die_x[i]) for i, C in enumerate(cmp)]
[Parameter('KD_Kx_A2_die_%s' % C.name, KD_Kx_A2_die[i]) for i, C in enumerate(cmp)]
[k_fate('k_A2_die_%s' % C.name, par['k_A2_die_0_%s' % C.name], par['k_A2_die_x_%s' % C.name],
        par['KD_Kx_A2_die_%s' % C.name], obs['Y_obs_%s' % C.name]) for i, C in enumerate(cmp)]
[Rule('A2_die_%s' % C.name, A2()**C >> None, exp['k_A2_die_%s' % C.name]) for i, C in enumerate(cmp)]

# #### Y div and death #####

[Parameter('k_Y_div_0_%s' % C.name, k_Y_div_0[i]) for i, C in enumerate(cmp)]
[Parameter('k_Y_div_x_%s' % C.name, k_Y_div_x[i]) for i, C in enumerate(cmp)]
[Parameter('KD_Kx_Y_div_%s' % C.name, KD_Kx_Y_div[i]) for i, C in enumerate(cmp)]
[k_fate('k_Y_div_%s' % C.name, par['k_Y_div_0_%s' % C.name], par['k_Y_div_x_%s' % C.name],
        par['KD_Kx_Y_div_%s' % C.name], obs['NE_all_%s' % C.name]) for i, C in enumerate(cmp)]
[Rule('Y_div_%s' % C.name, Y()**C >> Y()**C + Y()**C, exp['k_Y_div_%s' % C.name]) for i, C in enumerate(cmp)]

[Parameter('k_Y_die_%s' % C.name, k_Y_die[i]) for i, C in enumerate(cmp)]
[Rule('Y_die_%s' % C.name, Y()**C >> None, par['k_Y_die_%s' % C.name]) for i, C in enumerate(cmp)]

# Carrying capacity for total cell count
CC = 1e6

[Expression('k_A_cc_%s' % C.name, (exp['k_A_div_%s' % C.name] - exp['k_A_die_%s' % C.name])/CC)
 for i, C in enumerate(cmp)]
[Expression('k_N_cc_%s' % C.name, (exp['k_N_div_%s' % C.name] - exp['k_N_die_%s' % C.name])/CC)
 for i, C in enumerate(cmp)]
[Expression('k_A2_cc_%s' % C.name, (exp['k_A2_div_%s' % C.name] - exp['k_A2_die_%s' % C.name])/CC)
 for i, C in enumerate(cmp)]
[Expression('k_Y_cc_%s' % C.name, (exp['k_Y_div_%s' % C.name] - par['k_Y_die_%s' % C.name])/CC)
 for i, C in enumerate(cmp)]

# TODO: Report that in the Piecewise functions below, <= works but == doesn't
[Expression('rate_A_cc_%s' % C.name,
            Piecewise((0, obs['A_obs_%s' % C.name] <= 0),
                      (exp['k_A_cc_%s' % C.name]*Cells_TOT/obs['A_obs_%s' % C.name], True)))
    for i, C in enumerate(cmp)]
[Expression('rate_N_cc_%s' % C.name,
            Piecewise((0, obs['N_obs_%s' % C.name] <= 0),
                      (exp['k_N_cc_%s' % C.name]*Cells_TOT/obs['N_obs_%s' % C.name], True)))
    for i, C in enumerate(cmp)]
[Expression('rate_A2_cc_%s' % C.name,
            Piecewise((0, obs['A2_obs_%s' % C.name] <= 0),
                      (exp['k_A2_cc_%s' % C.name]*Cells_TOT/obs['A2_obs_%s' % C.name], True)))
    for i, C in enumerate(cmp)]
[Expression('rate_Y_cc_%s' % C.name,
            Piecewise((0, obs['Y_obs_%s' % C.name] <= 0),
                      (exp['k_Y_cc_%s' % C.name]*Cells_TOT/obs['Y_obs_%s' % C.name], True)))
    for i, C in enumerate(cmp)]

[Rule('A_cc_%s' % C.name, A()**C + A()**C >> A()**C, exp['rate_A_cc_%s' % C.name]) for C in cmp]
[Rule('N_cc_%s' % C.name, N()**C + N()**C >> N()**C, exp['rate_N_cc_%s' % C.name]) for C in cmp]
[Rule('A2_cc_%s' % C.name, A2()**C + A2()**C >> A2()**C, exp['rate_A2_cc_%s' % C.name]) for C in cmp]
[Rule('Y_cc_%s' % C.name, Y()**C + Y()**C >> Y()**C, exp['rate_Y_cc_%s' % C.name]) for C in cmp]

# #### Differentiation (state transitions) #####

[Parameter('kf_diff_A_N_%s' % C.name, kf_diff_A_N[i]) for i, C in enumerate(cmp)]
[Parameter('kr_diff_A_N_%s' % C.name, kr_diff_A_N[i]) for i, C in enumerate(cmp)]
[Rule('A_diff_N_%s' % C.name, A()**C | N()**C, par['kf_diff_A_N_%s' % C.name],
      par['kr_diff_A_N_%s' % C.name]) for i, C in enumerate(cmp)]

[Parameter('kf_diff_A_A2_%s' % C.name, kf_diff_A_A2[i]) for i, C in enumerate(cmp)]
[Parameter('kr_diff_A_A2_%s' % C.name, kr_diff_A_A2[i]) for i, C in enumerate(cmp)]
[Rule('A_diff_A2_%s' % C.name, A()**C | A2()**C, par['kf_diff_A_A2_%s' % C.name],
      par['kr_diff_A_A2_%s' % C.name]) for i, C in enumerate(cmp)]

[Parameter('kf_diff_N_A2_%s' % C.name, kf_diff_N_A2[i]) for i, C in enumerate(cmp)]
[Parameter('kr_diff_N_A2_%s' % C.name, kr_diff_N_A2[i]) for i, C in enumerate(cmp)]
[Rule('N_diff_A2_%s' % C.name, N()**C | A2()**C, par['kf_diff_N_A2_%s' % C.name],
      par['kr_diff_N_A2_%s' % C.name]) for i, C in enumerate(cmp)]

[Parameter('kf_diff_N_Y_%s' % C.name, kf_diff_N_Y[i]) for i, C in enumerate(cmp)]
[Parameter('kr_diff_N_Y_%s' % C.name, kr_diff_N_Y[i]) for i, C in enumerate(cmp)]
[Rule('N_diff_Y_%s' % C.name, N()**C | Y()**C, par['kf_diff_N_Y_%s' % C.name],
      par['kr_diff_N_Y_%s' % C.name]) for i, C in enumerate(cmp)]

# Invasion (epithelium to stroma)

Parameter('kf_A_epi_to_stroma', 0.1)  # 1e-2
Parameter('kf_N_epi_to_stroma', 0.1)  # 1e-2
Parameter('kf_A2_epi_to_stroma', 0.1)  # 1e-2
Parameter('kf_Y_epi_to_stroma', 0.1)  # 1e-2

Parameter('kr_A_epi_to_stroma', 0.02)  # 0
Parameter('kr_N_epi_to_stroma', 0.02)  # 1e-2
Parameter('kr_A2_epi_to_stroma', 0.02)  # 0
Parameter('kr_Y_epi_to_stroma', 0.02)  # 1e-2

Rule('A_epi_to_stroma', A()**E | A()**S, kf_A_epi_to_stroma, kr_A_epi_to_stroma)
Rule('N_epi_to_stroma', N()**E | N()**S, kf_N_epi_to_stroma, kr_N_epi_to_stroma)
Rule('A2_epi_to_stroma', A2()**E | A2()**S, kf_A2_epi_to_stroma, kr_A2_epi_to_stroma)
Rule('Y_epi_to_stroma', Y()**E | Y()**S, kf_Y_epi_to_stroma, kr_Y_epi_to_stroma)

# SIMULATIONS
sim = ScipyOdeSimulator(model, verbose=True)

# Simulate period before invasion begins
param_values = {'kf_A_epi_to_stroma': 0, 'kf_N_epi_to_stroma': 0,
                'kf_A2_epi_to_stroma': 0, 'kf_Y_epi_to_stroma': 0}
t_invasion = 100
tspan1 = np.linspace(0, t_invasion, t_invasion+1)
output1 = sim.run(tspan=tspan1, param_values=param_values)

# Simulation invasion
param_values = {'kf_A_epi_to_stroma': 0.1, 'kf_N_epi_to_stroma': 0.1,
                'kf_A2_epi_to_stroma': 0.1, 'kf_Y_epi_to_stroma': 0.1}
tspan2 = np.linspace(t_invasion, t_invasion+50, 51)
output2 = sim.run(tspan2, param_values=param_values)

fig, ax = plt.subplots(nrows=3, ncols=2, sharex='all', figsize=(9.6, 8.0))  # (6.4, 4.8)
t_left = 0.8*tspan1[-1]
t_right = tspan2[-1]
plt.xlim(left=t_left, right=t_right)

for lab, tspan, x in zip([True, False], [tspan1, tspan2], [output1, output2]):
# for lab, tspan, x in zip([True], [tspan2], [output2]):
    row = 0
    for C_name in [c.name for c in cmp] + ['TOT']:

        obs_name = ['A_obs_%s' % C_name,
                    'A2_obs_%s' % C_name,
                    'N_obs_%s' % C_name,
                    'Y_obs_%s' % C_name]

        color = ['darkred', 'r', 'c', 'b']

        col = 0  # cell counts
        label = []
        suffix = ''
        for i, name in enumerate(obs_name):
            label.append(name[:name.find('_')])
            if C_name == 'E':
                suffix = '_epith'
            elif C_name == 'S':
                suffix = '_stroma'
            else:
                suffix = '_total'
            label[-1] += suffix
            label_subtype = label[-1]
            if not lab:
                label_subtype = None
            ax[row, col].plot(tspan, x.all[name], lw=2, label=label_subtype, color=color[i])
        label_tot = 'Cells%s' % suffix
        if not lab:
            label_tot = None
        ax[row, col].plot(tspan, x.all['Cells_%s' % C_name], 'k--', lw=2, label=label_tot)
        if row == 2:
            ax[row, col].set_xlabel('time (d)')  # , fontsize=16)
        ax[row, col].set_ylabel('cell count')  # , fontsize=16)
        ax[row, col].xaxis.set_tick_params()  # labelsize=14)
        ax[row, col].yaxis.set_tick_params()  # labelsize=14)
        ax[row, col].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax[row, col].legend(loc=0, fontsize=6)

        col = 1  # ratios of cell counts
        cell_tot = sum(x.all[name] for name in obs_name)
        sum_prev = np.sum(np.array([val for val in x.all[obs_name[i]]]) for i in range(len(obs_name)-1))
        ax[row, col].fill_between(tspan, [1]*len(tspan), sum_prev / cell_tot, label=label[-1], color=color[-1])
        for i in range(len(obs_name) - 2, 0, -1):
            ax[row, col].fill_between(tspan, sum_prev / cell_tot, (sum_prev - x.all[obs_name[i]]) / cell_tot,
                                      label=label[i], color=color[i])
            sum_prev -= x.all[obs_name[i]]
        ax[row, col].fill_between(tspan, sum_prev / cell_tot, label=label[0], color=color[0])
        # RPM mice: 4, 5.7, 7.4, 9 weeks (RPM1-4, estimates)
        weeks_mice = np.linspace(4, 9, 4)
        for w in weeks_mice:
            ax[row, col].plot([t_invasion + (w-weeks_mice[0])*7+1, t_invasion + (w-weeks_mice[0])*7+1], [0, 1], 'w--', lw=2)
        ##########
        if row == 2:
            ax[row, col].set_xlabel('time (d)')  # , fontsize=16)
        ax[row, col].set_ylabel('cell fraction')  # , fontsize=16)
        ax[row, col].xaxis.set_tick_params()  # labelsize=14)
        ax[row, col].yaxis.set_tick_params()  # labelsize=14)
        # ax[row,col].legend(loc=(0.75, 0.6), framealpha=1)
        row += 1

plt.tight_layout()
plt.show()
