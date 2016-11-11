from __future__ import print_function, division
import numpy as np
from tabulate import tabulate

#T = 40.0
# = 60.0
#T = 80.0
#T = 100.0

T_0 = np.array([40.0, 60.0, 80.0, 100.0])

labels = ['Ideal gas', 'RK', 'SRK', 'PR', 'LK', 'SW']

mean_r = np.empty((len(T_0), len(labels)-1))
mean_c = np.empty((len(T_0), len(labels)-1))
max_r = np.empty((len(T_0), len(labels)-1))
max_c = np.empty((len(T_0), len(labels)-1))

for j, T in enumerate(T_0):
    rho = np.loadtxt('data/rho_{}.txt'.format(T))
    c = np.loadtxt('data/c_{}.txt'.format(T))

    for i, label in enumerate(labels[:-1]):
        d_r = 1 - rho[:, i+1]/rho[:, -1]
        d_c = 1 - c[:, i+1]/c[:, -1]

        mean_r[j, i] = np.mean(d_r)*100
        mean_c[j, i] = np.mean(d_c)*100

        max_r[j, i] = d_r[np.argmax(np.abs(d_r))]*100
        max_c[j, i] = d_c[np.argmax(np.abs(d_c))]*100

header_list = ['T'] + labels
header = ', '.join(header_list)

mean_r_tab = np.hstack([T_0[:, None], mean_r])
mean_c_tab = np.hstack([T_0[:, None], mean_c])
max_r_tab = np.hstack([T_0[:, None], max_r])
max_c_tab = np.hstack([T_0[:, None], max_c])

# Print tables
#######################
tablefmt = 'orgtbl'
'''
- "plain"
- "simple"
- "grid"
- "fancy_grid"
- "pipe"
- "orgtbl"
- "jira"
- "psql"
- "rst"
- "mediawiki"
- "moinmoin"
- "html"
- "latex"
- "latex_booktabs"
- "textile"
'''
#######################
print('Mean % difference in density')
print(tabulate(mean_r_tab, headers=header_list, tablefmt=tablefmt))
print()
print('Mean % difference in speed of sound')
print(tabulate(mean_c_tab, headers=header_list, tablefmt=tablefmt))
print()
print('Max. % difference in density')
print(tabulate(max_r_tab, headers=header_list, tablefmt=tablefmt))
print()
print('Max. % difference in speed of sound')
print(tabulate(max_c_tab, headers=header_list, tablefmt=tablefmt))

# save tables
np.savetxt('tables/rho_mean_diff.txt', mean_r_tab, header=header)
np.savetxt('tables/c_mean_diff.txt', mean_c_tab, header=header)
np.savetxt('tables/rho_max_diff.txt', max_r_tab, header=header)
np.savetxt('tables/c_max_diff.txt', max_c_tab, header=header)
