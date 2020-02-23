#!/usr/bin/python3
# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    sns.set(style="darkgrid") # set grid style of the plots

    numbers = [5, 10, 50, 100, 500, 1000, 2500, 5000, 7500, 10000] # length of sequences
    labels = ["measured"] * len(numbers) # set data label to measured
    # times_back = [2.5272369384765625e-05, 4.1961669921875e-05, 0.00023984909057617188, 0.0003578662872314453, 0.0017561912536621094, 0.003718852996826172, 0.009850740432739258, 0.04775285720825195, 0.048211097717285156, 0.0675346851348877]
    times_back = [2.6941299438476562e-05, 0.00013709068298339844, 0.0002129077911376953, 0.0003669261932373047, 0.0017430782318115234, 0.003720998764038086, 0.009377002716064453, 0.01937389373779297, 0.04137086868286133, 0.05624675750732422]

    for i in range(len(numbers ) + 1):
        times_back.append(times_back[0] / 5 * numbers[i])
        numbers.append(numbers[i])
        labels.append("expected")

    # plot the measurements of the trace back algorithm
    df1 = pd.DataFrame({'length of sequence': numbers, 'time [s]': times_back, 'label' : labels})
    plot1 = sns.lineplot(x = 'length of sequence', y = 'time [s]', hue = 'label', color = 'green', style = 'label', data = df1)
    plot1.set(ylim = (0, max(df1['time [s]'])))
    plt.show()

    numbers = [5, 10, 50, 100, 500, 1000, 2500, 5000, 7500, 10000] # length of sequences
    labels = ["measured"] * len(numbers) # set data label to measured
    # times_aff = [0.0003027915954589844, 0.0009229183197021484, 0.01734018325805664, 0.07471680641174316, 1.7508139610290527, 7.1654908657073975, 43.291638135910034, 186.4585838317871, 417.8935899734497, 712.7083909511566]
    times_aff = [0.0002791881561279297, 0.0012280941009521484, 0.016415834426879883, 0.06594514846801758, 1.6978788375854492, 6.503705024719238, 40.73264408111572, 166.2349147796631, 373.72168827056885, 677.711991071701]

    # add expected O(n^2) function and compare to measurements
    for i in range(len(numbers ) + 1):
        times_aff.append(times_aff[0] / 5 * numbers[i]**2 / 10)
        numbers.append(numbers[i])
        labels.append("expected")

    df2 = pd.DataFrame({'length of sequence': numbers, 'time [s]': times_aff, 'label' : labels})
    plot2 = sns.lineplot(x = 'length of sequence', y = 'time [s]', hue = 'label', style = 'label', data = df2)
    plot2.set(ylim = (0, max(df2['time [s]'])))
    plt.show()


if __name__ == '__main__':
    main()
