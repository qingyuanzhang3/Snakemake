import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

data = Table.read(snakemake.input[0], format="ascii")
time_data = np.array(data["year"])
d14c_data = np.array(data["d14c"])
d14c_data_error = np.array(data["sig_d14c"])
plt.errorbar(time_data, d14c_data, yerr=d14c_data_error, fmt='.k')
plt.ylabel("$\Delta^{14}$C (‰)")
plt.xlabel("Calender years")
plt.title(snakemake.input[0])
plt.savefig(snakemake.output[0])
