import matplotlib.pyplot as plt

def customize_matplotlib():
    font_size = 14

    rc_params = { 
                'font.size': font_size,
                #'font.family': 'Times New Roman',
                'axes.labelsize': font_size,
                'axes.titlesize': 1.2*font_size,
                'axes.grid': True, 
                'xtick.labelsize': font_size,
                'ytick.labelsize': font_size,
                #'text.usetex': True,
                'legend.fontsize': font_size,
                #'legend.frameon': False
                }
    plt.rcParams.update(rc_params)
