import argparse
import random
import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from eopr import get_ubar, error_bars

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Ellipsoidal Optimal Recovery')
    parser.add_argument('--data', type=str, default='../data/basque_abadi_reshaped.csv', help='data file name')
    parser.add_argument('--pretreatment_end', type=str, default="1971", help='index int of column')
    parser.add_argument('--results', type=str, default='results', help='results output directory')
    parser.add_argument('--figures', type=str, default='figures', help='figures output directory')

    args = parser.parse_args()

    datafile = args.data
    pretreatment_end = args.pretreatment_end
    results_dir = args.results
    figures_dir = args.figures

    if not os.path.exists(figures_dir):
        os.mkdir(figures_dir)
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    
    # first row is the treated row
    df = pd.read_csv(datafile, index_col=0)
    cols = np.array(df.columns)
    TrainingEnd = np.where(cols==pretreatment_end)[0][0] # index
    M = np.array(df)
    m = M[0] # treated unit
    M = M[1:] # control units 
    

    vals = []
    op_errors = []
    regs = np.arange(0.01,1, 0.01)
    xx = np.array([i for i in range(TrainingEnd)])
    for reg in regs:
        xx,ubar,PHI, R, Q = get_ubar(M,m,xx,reg)
        vals.append(np.sqrt(np.mean((ubar[:TrainingEnd] - m[:TrainingEnd])**2)))

    mm = np.argmin(vals)
    min_reg = regs[mm]
    xx,ubar,PHI, R, Q = get_ubar(M,m,xx,min_reg)
    
    # save estimates
    pd.DataFrame({'estimates':ubar}).to_csv(f'{results_dir}/estimates.csv')

    # get error bars
    max_error = error_bars(xx,ubar,M, PHI, R, Q)

    # visualize
    fig, ax = plt.subplots(figsize=(8,6))
    lw = 2.5
    fs = 20
    plt.plot(m, '-',label='treated unit',fillstyle='none',color='black',linewidth=lw-0.5)

    plt.plot(cols, ubar, '--',label='EOpR',fillstyle='none',color='tab:blue',linewidth=lw,zorder=10)
    plt.axvline(str(int(pretreatment_end)-1), linestyle = '--', color = 'black',linewidth=1,label='intervention')
    plt.xticks(rotation=90,fontsize=fs-8)
    plt.yticks(fontsize=fs-8)

    plt.ylabel('outcome',fontsize=fs-8)
    plt.xlabel('Year',fontsize=fs-5)
    plt.legend(fontsize=fs-9)
    plt.tight_layout()
    plt.savefig(f'{figures_dir}/estimation_viz.png',dpi=300)

    # visualize error bars
    fig, ax = plt.subplots(figsize=(8,6))
    lw = 2.5
    fs = 20
    plt.plot(m, '-',label='treated unit',fillstyle='none',color='black',linewidth=lw-0.5)
    
    plt.plot(cols, ubar, '--',label='EOpR',fillstyle='none',color='tab:blue',linewidth=lw,zorder=10)
    plt.plot((ubar.T + max_error)[0], '--',label='worst case',fillstyle='none',linewidth=lw-1,color='gray')
    plt.plot((ubar.T - max_error)[0], '--',fillstyle='none',linewidth=lw-1,color='gray')

    plt.axvline(str(int(pretreatment_end)-1), linestyle = '--', color = 'black',linewidth=1,label='intervention')

    plt.xticks(rotation=90,fontsize=fs-8)
    plt.yticks(fontsize=fs-8)

    plt.ylabel('outcome',fontsize=fs-8)
    plt.xlabel('Year',fontsize=fs-5)
    plt.legend(fontsize=fs-9)
    plt.tight_layout()
    plt.savefig(f'{figures_dir}/estimation_with_errors_viz.png',dpi=300)

