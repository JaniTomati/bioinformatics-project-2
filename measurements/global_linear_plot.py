#!/usr/bin/python3
# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    sns.set(style="darkgrid") # set grid style of the plots

    numbers = [5, 10, 50, 100, 500, 1000, 2500, 5000, 7500, 10000] # length of sequences
    labels = ["measured"] * len(numbers) # set data label to measured

    # plot the measurements of the trace back algorithm
    times_back = [2.5033950805664062e-05, 9.608268737792969e-05, 0.00016498565673828125, 0.00032806396484375, 0.0016469955444335938, 0.0032930374145507812, 0.008365154266357422, 0.016770124435424805, 0.02606797218322754, 0.034770965576171875]
    df1 = pd.DataFrame({'length of sequence': numbers, 'time [s]': times_back, 'label' : labels})
    plot1 = sns.lineplot(x = 'length of sequence', y = 'time [s]', hue = 'label', color = 'green', data = df1)
    plot1.set(ylim = (0, max(df1['time [s]'])))
    plt.show()

    # add expected O(n^2) function and compare to measurements
    times_opt = [0.00017786026000976562, 0.0010161399841308594, 0.011538982391357422, 0.03510594367980957, 0.8658099174499512, 3.5476551055908203, 22.301820039749146, 91.95854210853577, 204.99176812171936, 360.73924803733826]

    for i in range(len(numbers ) + 1):
        times_opt.append(times_opt[0] / 5 * numbers[i]**2 / 10)
        numbers.append(numbers[i])
        labels.append("expected")

    df2 = pd.DataFrame({'length of sequence': numbers, 'time [s]': times_opt, 'label' : labels})
    plot2 = sns.lineplot(x = 'length of sequence', y = 'time [s]', hue = 'label', style = 'label', data = df2)
    plot2.set(ylim = (0, max(df2['time [s]'])))
    plt.show()


if __name__ == '__main__':
    main()
