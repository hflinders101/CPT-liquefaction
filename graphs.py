import pandas as pd
import matplotlib.pyplot as plt
import math
import os, glob
from functions import soil_parameters


def ordMag(number):
    if number == 0:
        return 1
    elif int(str(abs(number))[:1]) == 9:
        return math.floor(math.log(abs(number), 10)) + 1
    else:
        return math.floor(math.log(abs(number), 10))


def plotGridlines(max, min):
    small_rounding_dict = {
        -13: .00000000000001,
        -12: .0000000000001,
        -11: .000000000001,
        -10: .00000000001,
        -9: .0000000001,
        -8: .000000001,
        -7: .00000001,
        -6: .0000001,
        -5: .000001,
        -4: .00001,
        -3: .0001,
        -2: .001,
        -1: .01,
        0: .1,
        1: 1,
        2: 10,
        3: 100,
        4: 1000,
        5: 10000
    }
    large_rounding_dict = {
        -13: .1,
        -12: .1,
        -11: .1,
        -10: .1,
        -9: .1,
        -8: .1,
        -7: .1,
        -6: .1,
        -5: .1,
        -4: .1,
        -3: .1,
        -2: .1,
        -1: .1,
        0: 1,
        1: 1,
        2: 10,
        3: 100,
        4: 1000,
        5: 1000,
    }

    if (ordMag(max) < 0 or ordMag(min) < 0) and (ordMag(abs(max)) >= -0 or ordMag(abs(min)) >= -0):
        roundDown = large_rounding_dict.get(ordMag(min))
        roundUp = small_rounding_dict.get(ordMag(max))

    elif ordMag(max) < 0 or ordMag(min) < 0:
        roundDown = small_rounding_dict.get(ordMag(min))
        roundUp = small_rounding_dict.get(ordMag(max))
    else:
        roundDown = large_rounding_dict.get(ordMag(min))
        roundUp = large_rounding_dict.get(ordMag(max))

    upper_bound = ((max / roundUp + 9) // 10 * 10) * roundUp
    lower_bound = ((min / roundDown) // 10 * 10) * roundDown

    difference = upper_bound - lower_bound

    prime_numbers = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

    for prime in prime_numbers:
        if difference % prime == 0 or 5:
            difference = difference / prime
            break
    return [lower_bound, lower_bound + (difference) / 2, lower_bound + difference,
            lower_bound + difference + difference / 2, lower_bound + difference * 2, upper_bound]


def makePlots(depth_column, dependent_variables_list, counter=1):
    plot_name = os.path.basename(filename)

    columns = dependent_variables_list
    plots = dependent_variables_list

    if len(plots) == 1:
        fig, ax = plt.subplots()
        y = depth_column
        x = dependent_variables_list[0]
        plt.plot(x, y)
        plt.xlabel(dependent_variables_list[0].name)
        plt.ylabel(depth_column.name)
        plt.title(plot_name)
        fig.set_figheight(6.9)
        fig.set_figwidth(2.4)
        min = dependent_variables_list[0].min()
        if min == -9999:
            min = dependent_variables_list[0].nsmallest(2).iloc[-1]
        max = dependent_variables_list[0].max()
        ax.set_xticks(plotGridlines(max, min))
        ax.set_xlim(plotGridlines(max, min)[0], plotGridlines(max, min)[-1])
        ax.grid()
        ax.set_ylim(ymin=0)  # sets the start and end of graph to 0 and last depth
        ax.invert_yaxis()  # inverts y axis
        ax.tick_params(axis='x', rotation=90)
        fig.set_tight_layout(True)
        fig.savefig(filename + "\\" + plot_name + "_" + str(counter), bbox_inches='tight', dpi=300)
        plt.close(fig)

    else:
        fig1, (plots) = plt.subplots(1, len(dependent_variables_list), sharey=True)
        fig1.suptitle(plot_name)
        fig1.set_figheight(6.9)
        fig1.set_figwidth(12 / (5 / len(plots)))
        fig1.text(0.075 / (5 / (len(plots) * .5)), 0.5, depth_column.name, va='center', rotation='vertical')

        for column, num_plot in zip(columns, plots):

            num_plot.plot(column, depth_column)

            min = column.min()
            if min == -9999:
                min = column.nsmallest(2).iloc[-1]
            max = column.max()

            num_plot.set_xticks(plotGridlines(max, min))
            num_plot.set_xlim(plotGridlines(max,min)[0],plotGridlines(max,min)[-1])
            num_plot.set_title(column.name)

        for plot in plots:
            plot.tick_params(axis="x", rotation=90)
            plot.grid()

        plots[0].set_ylim(ymin=0)  # sets the start and end of graph to 0 and last depth
        plots[0].invert_yaxis()  # inverts y axis
        figure1 = fig1

        figure1.savefig(filename + "\\" + plot_name + "_" + str(counter), bbox_inches='tight', dpi=300)
        plt.close(figure1)


folder_path = r"C:\Users\jdundas2\Documents\paper tests"
for filename in glob.glob(os.path.join(folder_path, "*.xls*")):
    print(filename)
    df = pd.read_excel(filename)
    df = soil_parameters(df)
    # filename = filename.replace("cliquest homo", "cliquest homo graphed")
    filename = filename.removesuffix(".xls")
    filename = filename.removesuffix(".xlsx")
    plot_name = os.path.basename(filename)
    os.mkdir(filename)

    # Input the number of figures and graphs per figure to create
    num_of_figures = 6
    graphs_per_figure = 5

    # Use this line to select a range of columns
    variables = df.columns.values.tolist()[1:30]

    # Use these two lines to select individual columns by index number
    # df_custom = df.iloc[:, [3, 5, 1, 7]]
    # variables = df_custom.columns.values.tolist()

    # Use this line to select individual columns by name
    # variables = ['u (kPa)']

    total_num_variables = len(variables)
    count = 0

    inputs = []

    for x in range(num_of_figures):
        inputs.append('input' + str(x))

    for i in range(len(inputs)):
        if total_num_variables - count <= graphs_per_figure:
            remainder = total_num_variables - count
            dataframe = df[variables[count:count + remainder]]
        else:
            dataframe = df[variables[count:count + graphs_per_figure]]

        inputs[i] = []
        for col in dataframe:
            inputs[i].append(df[col])
        count += graphs_per_figure

    for i in range(num_of_figures):
        counter = i + 1
        makePlots(df['Depth (m)'], inputs[i], counter)
