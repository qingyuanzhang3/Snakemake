import matplotlib.pyplot as plt
import numpy as np
import ticktack
from ticktack import fitting
from matplotlib.lines import Line2D
import matplotlib as mpl
import seaborn as sns
mpl.style.use('seaborn-colorblind')


events = ["775AD-Sharp", "775AD-Prolonged", "993AD", "663BCE", "5259BCE", "5410BCE", "7176BCE"]
titles = ["775CE Sharp-Rise", "775CE Prolonged-Rise", "993CE", "663BCE", "5259BCE", "5410BCE", "7176BCE"]
cbm_models = ["Guttler15", "Buntgen18", "Brehm21",]
steady_state = []
for cbm_model in cbm_models:
    cbm = ticktack.load_presaved_model(cbm_model, production_rate_units='atoms/cm^2/s')
    sf = fitting.SingleFitter(cbm, cbm_model=cbm_model, hemisphere="north")
    steady_state.append(sf.steady_state_production)

cbm_names = ["Güttler et al, 2015: 11-box", "Büntgen et al, 2018: 22-box", "Brehm et al, 2021: 22-box"]
colors = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#e41a1c', '#56B4E9', '#999999', ]

Guttler15 = np.zeros((12000, len(events)))
Buntgen18 = np.zeros((12000, len(events)))
Brehm21 = np.zeros((12000, len(events)))

units = "atoms"
if units == "atoms":
    conversion = 1
else:
    conversion = 14.003242 / 6.022 * 5.11 * 31536. / 1e5
for i, event in enumerate(events):
    Guttler15[:, i] = np.load("chain/{}_Guttler15.npy".format(event))[:, 4]
    Guttler15[:, i] = 10**Guttler15[:, i] * conversion
for i, event in enumerate(events):
    Buntgen18[:, i] = np.load("chain/{}_Buntgen18.npy".format(event))[:, 4]
    Buntgen18[:, i] = 10 ** Buntgen18[:, i] * conversion
for i, event in enumerate(events):
    Brehm21[:, i] = np.load("chain/{}_Brehm21.npy".format(event))[:, 4]
    Brehm21[:, i] = 10 ** Brehm21[:, i] * conversion

fmt = ["-", "--", ":"]
custom_lines = [Line2D([0], [0], ls=fmt[i], color="k", lw=1.5, label=cbm_names[i]) for i in range(len(cbm_names))]
fig = plt.figure(figsize=(10, 5), dpi=200, constrained_layout=True)
spec = fig.add_gridspec(ncols=1, nrows=1)
ax = fig.add_subplot(spec[0, 0])
for i in range(0, 7):
    sns.kdeplot(Guttler15[:, i] / steady_state[0], ls=fmt[0], color=colors[i], ax=ax, clip=(0, 15));
    sns.kdeplot(Buntgen18[:, i] / steady_state[1], color=colors[i], ls=fmt[1], ax=ax, clip=(0, 15));
    sns.kdeplot(Brehm21[:, i] / steady_state[2], ls=fmt[2], color=colors[i], ax=ax, clip=(0, 15));
ax.legend(handles=custom_lines, frameon=False, fontsize=10, loc="upper left");

plt.text(x=2.65, y=4, s=titles[0], color=colors[0])
plt.text(x=2.65, y=5, s=titles[1], color=colors[1])
plt.text(x=1.3, y=2.45, s=titles[2], color=colors[2])
plt.text(x=2.2, y=5, s=titles[3], color=colors[3])
plt.text(x=4.2, y=3, s=titles[4], color=colors[4])
plt.text(x=0.7, y=3.6, s=titles[5], color=colors[5])
plt.text(x=3.7, y=4.7, s=titles[6], color=colors[6])


fig.supxlabel("spike production (years of steady state production)", fontsize=16, fontfamily="serif", fontweight="roman");
fig.savefig(snakemake.output[0])
