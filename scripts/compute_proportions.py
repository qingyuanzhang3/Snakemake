import numpy as np
from ticktack import fitting

for path in snakemake.input:
    print(snakemake.input)
    chain = np.load(path)
    p1 = chain[:, 2][10 ** chain[:, 2] < 1 / 2].size / chain[:, 2].size
    p2 = chain[:, 2][(10 ** chain[:, 2] > 1 / 2) & (10 ** chain[:, 2] < 1)].size / chain[:, 2].size
    p3 = chain[:, 2][10 ** chain[:, 2] > 1].size / chain[:, 2].size

    with open('proportions.txt', 'a') as f:
        f.write(path + "\n")
        f.write("<0.5 | 0.5 - 1 | > 1 : %.2f | %.2f | %.2f \n" % (p1, p2, p3))
        f.write("\n")
