import matplotlib.pyplot as plt
import numpy as np
import ticktack
from ticktack import fitting
from matplotlib.lines import Line2D
import matplotlib as mpl
import os
mpl.style.use('seaborn-colorblind')


cbm_models = ["Guttler15", "Buntgen18", "Brehm21",]
cbm_names = ["Güttler et al, 2015: 11-box", "Büntgen et al, 2018: 22-box", "Brehm et al, 2021: 22-box"]
colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
custom_lines = [Line2D([0], [0], color=colors[i], lw=1.5, label=cbm_names[i]) for i in range(len(cbm_models))]
custom_lines.append(plt.errorbar([0], [0], yerr=1, fmt="ok", capsize=3, label="average $\Delta^{14}$C"))
custom_lines.append(Line2D([0], [0], ls="--", color="k", marker="o", lw=1.5, label="average $\Delta^{14}$C NH"))
custom_lines.append(Line2D([0], [0], ls="-", color="k", marker="s", lw=1.5, label="average $\Delta^{14}$C SH"))
events = ["663BCE", "5259BCE", "775AD-Sharp", "5410BCE", "7176BCE"]
titles = ["663BCE", "5259BCE", "775CE Sharp-Rise", "5410BCE", "7176BCE"]
intervals = [None, 5, None, 5, 5]

def remove_frame_top(ax):
    ax.get_xaxis().set_visible(False);
    ax.spines['bottom'].set_visible(False)

def remove_frame_bot(ax):
    ax.spines['top'].set_visible(False)

def remove_frame(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.xaxis.set_ticks([]); ax.yaxis.set_ticks([])
    ax.xaxis.set_ticklabels([]); ax.yaxis.set_ticklabels([])

def plot_CP(ax0, ax1, event=None, title=None, interval=None, markersize=6, markersize2=3, capsize=3, elinewidth=3):
    soln_path = ["non-parametric/solutions/{}_{}.npy".format(event, cbm_model) for cbm_model in cbm_models]
    chain_path = ["non-parametric/chain/{}_{}.npy".format(event, cbm_model) for cbm_model in cbm_models]
    inverse_solver_path = ["non-parametric/solver/{}_{}.npy".format(event, cbm_model) for cbm_model in cbm_models]
    chain = np.load(inverse_solver_path[0])
    m, n = chain.shape
    merged_chain = np.zeros((m * 3, n))
    for i, IS in enumerate(inverse_solver_path):
        chain = np.load(inverse_solver_path[i])
        chain = np.where(chain < 0, 0, chain)
        merged_chain[i*m:(i+1)*m, :] = chain
    fitting.plot_ControlPoints(average_path="data-CP/means/{}.csv".format(event),
                               soln_path=soln_path,
                               chain_path=chain_path,
                               merged_inverse_solver=merged_chain,
                               cbm_models=cbm_models,
                               hemisphere="north",
                               directory_path="data-CP/{}".format(event),
                               labels=False, axs=[ax0, ax1], interval=interval, markersize=markersize,
                               markersize2=markersize2, capsize=capsize, elinewidth=elinewidth)
    ax0.set_title(title, fontsize=14, fontfamily="serif",fontweight="roman")

def sp_event_775AD(ax0, ax1, title=None, capsize=6, markersize=3, markersize2=3, elinewidth=3):
    events = ["775AD-Prolonged-N", "775AD-Prolonged-S"]
    hemisphere = ["north", "south"]
    colors_prod = ["b", "r"]
    fmts = [["ok", "--k"], ["sk", "-k"]]
    size = 100
    size2 = 30

    for k, event in enumerate(events):
        inverse_solver_path = ["non-parametric/solver/{}_{}.npy".format(event, cbm_model) for cbm_model in cbm_models]
        chain = np.load(inverse_solver_path[0])
        m, n = chain.shape
        merged_chain = np.zeros((m * 3, n))
        for i, IS in enumerate(inverse_solver_path):
            chain = np.load(inverse_solver_path[i])
            chain = np.where(chain < 0, 0, chain)
            merged_chain[i*m:(i+1)*m, :] = chain
        sf = fitting.SingleFitter("Guttler15", cbm_model="Guttler15", hemisphere=hemisphere[k])
        sf.load_data("data-CP/means/{}.csv".format(event))
        ax1.errorbar(sf.time_data + sf.time_offset, np.median(chain, axis=0), fmt=colors_prod[k], drawstyle="steps", alpha=0.2)
        ax1.fill_between(sf.time_data + sf.time_offset, np.percentile(chain, 32, axis=0),
                         np.percentile(chain, 68, axis=0), step='pre', alpha=0.1,
                             color=colors_prod[k], edgecolor="none",
                             lw=1.5)

    for i, event in enumerate(events):
        for j, model in enumerate(models):
            cbm = ticktack.load_presaved_model(model, production_rate_units = 'atoms/cm^2/s')
            sf = fitting.SingleFitter(cbm, cbm_model=model, hemisphere=hemisphere[i])
            sf.load_data("data-CP/means/" + event + ".csv")
            sf.compile_production_model(model="control_points")
            soln = np.load("non-parametric/solutions/775AD-Prolonged_" + model + ".npy", allow_pickle=True)
            chain = np.load("non-parametric/chain/775AD-Prolonged_" + model + ".npy", allow_pickle=True)
            mu = np.mean(chain, axis=0)

            ax0.plot(sf.time_data_fine, sf.dc14_fine(mu), color=colors[j])
            ax1.plot(sf.control_points_time_fine, sf.interp_gp(sf.control_points_time_fine, mu), color=colors[j])
            idx = np.random.randint(len(chain), size=30)
            for param in chain[idx]:
                ax1.plot(sf.control_points_time_fine, sf.interp_gp(sf.control_points_time_fine, param),
                         alpha=0.2, color=colors[j])

        ax0.errorbar(sf.time_data + sf.time_offset, sf.d14c_data, yerr=sf.d14c_data_error, fmt=fmts[i][0], capsize=capsize,
                     markersize=markersize, elinewidth=elinewidth, label="average $\Delta^{14}$C", alpha=1)
        ax0.plot(sf.time_data + sf.time_offset, sf.d14c_data, fmts[i][1], alpha=1)
        file_names = [f for f in os.listdir("data-CP/" + event) if os.path.isfile(os.path.join("data-CP/" + event, f))]
        for file in file_names:
            sf.load_data("data-CP/" + event + '/' + file)
            ax0.errorbar(sf.time_data + sf.time_offset, sf.d14c_data, fmt="o", color="gray", yerr=sf.d14c_data_error, capsize=3,
                         alpha=0.2)
    ax0.set_title(title, fontsize=14, fontfamily="serif",fontweight="roman")
    ax1.set_xlim(sf.start, sf.end+1.2);

def sp_event_993AD(ax0, ax1, title=None, sp_event=None, capsize=6, markersize=3, markersize2=3, elinewidth=3):
    events = ["993AD-N", "993AD-S"]
    hemisphere = ["north", "south"]
    colors_prod = ["b", "r"]
    fmts = [["ok", "--k"], ["sk", "-k"]]
    size = 100
    size2 = 30

    for k, event in enumerate(events):
        inverse_solver_path = ["non-parametric/solver/{}_{}.npy".format(event, cbm_model) for cbm_model in cbm_models]
        chain = np.load(inverse_solver_path[0])
        m, n = chain.shape
        merged_chain = np.zeros((m * 3, n))
        for i, IS in enumerate(inverse_solver_path):
            chain = np.load(inverse_solver_path[i])
            chain = np.where(chain < 0, 0, chain)
            merged_chain[i*m:(i+1)*m, :] = chain
        sf = fitting.SingleFitter("Guttler15", cbm_model="Guttler15", hemisphere=hemisphere[k])
        sf.load_data("data-CP/means/{}.csv".format(event))
        ax1.errorbar(sf.time_data + sf.time_offset, np.median(chain, axis=0), fmt=colors_prod[k], drawstyle="steps", alpha=0.2)
        ax1.fill_between(sf.time_data + sf.time_offset, np.percentile(chain, 32, axis=0),
                         np.percentile(chain, 68, axis=0), step='pre', alpha=0.1,
                             color=colors_prod[k], edgecolor="none",
                             lw=1.5)

    for i, event in enumerate(events):
        for j, model in enumerate(models):
            cbm = ticktack.load_presaved_model(model, production_rate_units = 'atoms/cm^2/s')
            sf = fitting.SingleFitter(cbm, cbm_model=model, hemisphere=hemisphere[i])
            sf.load_data("data-CP/means/" + event + ".csv")
            sf.compile_production_model(model="control_points")
            soln = np.load("non-parametric/solutions/993AD_" + model + ".npy", allow_pickle=True)
            chain = np.load("non-parametric/chain/993AD_" + model + ".npy", allow_pickle=True)

            if event == "993AD-N":
                time_data = sf.time_data
                time_data_fine = sf.time_data_fine
                control_points_time = sf.control_points_time
                annual = sf.annual
                start = sf.start
                end = sf.end
            if event == "993AD-S":
                sf.time_data = time_data
                sf.time_data_fine = time_data_fine
                sf.control_points_time = control_points_time
                sf.annual = annual

            mu = np.mean(chain, axis=0)

            ax0.plot(sf.time_data_fine, sf.dc14_fine(mu), color=colors[j])
            ax1.plot(sf.control_points_time_fine, sf.interp_gp(sf.control_points_time_fine, mu), color=colors[j])
            idx = np.random.randint(len(chain), size=30)
            for param in chain[idx]:
                ax1.plot(sf.control_points_time_fine, sf.interp_gp(sf.control_points_time_fine, param),
                         alpha=0.2, color=colors[j])

        sf.load_data("data-CP/means/" + event + ".csv")
        ax0.errorbar(sf.time_data + sf.time_offset, sf.d14c_data, yerr=sf.d14c_data_error, fmt=fmts[i][0], capsize=capsize,
                     markersize=markersize, elinewidth=elinewidth, label="average $\Delta^{14}$C", alpha=1)
        ax0.plot(sf.time_data + sf.time_offset, sf.d14c_data, fmts[i][1], alpha=1)
        file_names = [f for f in os.listdir("data-CP/" + event) if os.path.isfile(os.path.join("data-CP/" + event, f))]
        for file in file_names:
            sf.load_data("data-CP/" + event + '/' + file)
            ax0.errorbar(sf.time_data + sf.time_offset, sf.d14c_data, fmt="o", color="gray", yerr=sf.d14c_data_error, capsize=3,
                         alpha=0.2)


    ax0.set_title(title, fontsize=14, fontfamily="serif",fontweight="roman")
    ax1.set_xlim(start, end+1);

models = ["Guttler15", "Brehm21", "Buntgen18"]
fig = plt.figure(figsize=(12, 24), dpi=200, constrained_layout=True)
fig.subplots_adjust(hspace=0.15)
spec = fig.add_gridspec(ncols=2, nrows=15, height_ratios=[1, 1, 1, 0.5] * 3 + [1, 1, 1])

# left 1
ax0 = fig.add_subplot(spec[0:2, 0])
ax1 = fig.add_subplot(spec[2, 0], sharex=ax0)
ax0.xaxis.set_visible(False);

# right 1
ax2 = fig.add_subplot(spec[0:2, 1])
ax3 = fig.add_subplot(spec[2, 1], sharex=ax2)
ax2.xaxis.set_visible(False);

# invisible block
fig.add_subplot(spec[3, :]).set_visible(False)

# second row
ax4 = fig.add_subplot(spec[4:6, 0])
ax5 = fig.add_subplot(spec[6, 0], sharex=ax4)
ax4.xaxis.set_visible(False);

ax6 = fig.add_subplot(spec[4:6, 1])
ax7 = fig.add_subplot(spec[6, 1], sharex=ax6)
ax6.xaxis.set_visible(False);

# invisible
fig.add_subplot(spec[7, :]).set_visible(False)

# third row
ax8 = fig.add_subplot(spec[8:10, 0])
ax9 = fig.add_subplot(spec[10, 0], sharex=ax8)
ax8.xaxis.set_visible(False);

ax10 = fig.add_subplot(spec[8:10, 1])
ax11 = fig.add_subplot(spec[10, 1], sharex=ax10)
ax10.xaxis.set_visible(False);

# invisible
fig.add_subplot(spec[11, :]).set_visible(False)

# fourth row
ax12 = fig.add_subplot(spec[12:14, 0])
ax13 = fig.add_subplot(spec[14, 0], sharex=ax12)
ax12.xaxis.set_visible(False);

ax14 = fig.add_subplot(spec[12:14, 1])
ax15 = fig.add_subplot(spec[14, 1], sharex=ax14)
ax14.legend(handles=custom_lines, frameon=False, loc="center", fontsize=14);
ax14.axis("off")
remove_frame(ax15)

# sup
ax13.set_xlabel("Year (BCE)", fontsize=16, fontfamily="serif", fontweight="roman")
ax15.set_xlabel("Year (CE)", fontsize=16, fontfamily="serif", fontweight="roman")
fig.supylabel("Production rate (atoms/cm$^2$/s) | $\Delta^{14}$C (‰)", x=0.05, y=0.5,
              fontsize=16, fontfamily="serif", fontweight="roman");

# add plots
axs = [(ax0, ax1), (ax4, ax5), (ax6, ax7), (ax8, ax9), (ax12, ax13)]
sp_event_775AD(ax2, ax3, title="775CE Prolonged-Rise", markersize=5, markersize2=4, capsize=3, elinewidth=2)
sp_event_993AD(ax10, ax11, title="993CE", markersize=5, markersize2=4, capsize=3, elinewidth=2)
for i, ax in enumerate(axs):
    if events[i] == "7176BCE":
        plot_CP(ax[0], ax[1], event=events[i], title=titles[i], interval=intervals[i], markersize=3,
                    markersize2=3, capsize=2, elinewidth=2)
    else:
        plot_CP(ax[0], ax[1], event=events[i], title=titles[i], interval=intervals[i], markersize=5,
                markersize2=4, capsize=3, elinewidth=2)
fig.savefig(snakemake.output[0], bbox_inches='tight')
