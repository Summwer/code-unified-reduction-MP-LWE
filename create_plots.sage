##########################################
## This file contains functions         ##
## to create the figures displayed      ##
## in the article                       ##
## from the data obtained after running ##
## the tests from run_tests.sage        ##
##########################################

## Note: loading this file will not directly create the plots, you should uncomment the corresponding lines below to create the plots

## Requirement: when running this code, you should execute the code in a repository containing a sub-repository named "figure" (the figures obtained will be saved there)

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import colormaps
import numpy as np

load("run_tests.sage")

########################################
## Some auxilliary plotting functions ##
########################################

## This function takes a input a list of data (e.g., A or B), an offset (which measures
## how the boxplot is moved on the x axis), the width of the boxplots, and colors for
## the edges and for the inside of the box
def draw_plot(data, offset, width, edge_color, fill_color):
    pos = np.arange(len(data))+offset
    bp = plt.boxplot(data, positions= pos, widths=width, patch_artist=True, manage_ticks=False, showfliers = False, whis = (10,90))
    for element in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)
    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)
        
## list of colors for pictures
list_colors = colormaps['tab20'].colors
dark_colors = [list_colors[i] for i in range(len(list_colors)) if i%2 == 0]
light_colors = [list_colors[i] for i in range(len(list_colors)) if i%2 == 1]

## list of possible markers
list_markers = ['o', 'v','^','s','p','*']

def format_B_label(B):
    #if B > 0:
        #expon = int(np.round(np.log2(B)))
        #if 2**expon == B:
        #    return r'$2^{%d}$' % expon
    return str(B)

#################################
###  Plots for D(f) quantity  ###
#################################

## Figure for fixed B and increasing degree

def plot_Df_unif_varying_degree():
    plt.close() ## reset plot

    ### loading data for B = 10
    data_10 = []
    x_labels = []
    print("loading data for D(f) with varying degree...")
    for (degree,B,nb_tests,_,__) in list_input_Df_unif:
        if B == 10:
            x_labels += [degree]
            load("data/Df_"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+"_unif.sage")
            data_10 += [[x[0] for x in res]]
            
    ### loading data for B = 1000
    data_1000 = []
    for deg in x_labels:
        for (degree,B,nb_tests,_,__) in list_input_Df_unif:
            if B == 1000 and degree == deg:
                load("data/Df_"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+"_unif.sage")
                data_1000 += [[x[0] for x in res]]
                break
    assert(len(data_10) == len(data_1000))
    print("....done")

    ## create the axes
    fig, ax = plt.subplots()

    ## set the name of the axes
    ax.set_ylabel('D(f)')
    ax.set_xlabel('Degree of f')
    
    ## plot a legend
    plt.legend(handles = [
        Patch(facecolor = light_colors[0], edgecolor=dark_colors[0], label = "B = 10"),
        Patch(facecolor = light_colors[1], edgecolor=dark_colors[1], label = "B = 1000"),
    ], loc = 'upper left')

    ## position the labels appearing on the x axis
    plt.xticks(range(len(x_labels)),x_labels)

    ## plot the boxplots for A and B with differenc colors and offset
    draw_plot(data_10, -0.15, 0.25, dark_colors[0], light_colors[0])
    draw_plot(data_1000, +0.15,0.25, dark_colors[1], light_colors[1])

    plt.savefig("figures/Df_unif_varying_degree.pdf")
 
    
## Figure for fixed (smallish) degree and varying B
def plot_Df_unif_varying_B():
    plt.close() ## reset plot

    ### loading data for degree 50
    data_50 = []
    data_100 = []
    x_labels = [10^i for i in range(1,9)]
    print("loading data for D(f) with varying bound B...")
    for B in x_labels:
        (degree,nb_tests) = (50,1000)
        load("data/Df_"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+"_unif.sage")
        data_50 += [[x[0] for x in res]]
        
        (degree,nb_tests) = (100,1000)
        load("data/Df_"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+"_unif.sage")
        data_100 += [[x[0] for x in res]]
    print("....done")

    ## create the plot
    fig, ax = plt.subplots()

    ## set the name of the axes
    ax.set_ylabel('D(f)')
    ax.set_xlabel('Bound B on the coefficients')
    
    ## plot a legend
    plt.legend(handles = [
        Patch(facecolor = light_colors[2], edgecolor=dark_colors[2], label = "degree 50"),
        Patch(facecolor = light_colors[3], edgecolor=dark_colors[3], label = "degree 100"),
    ], loc = 'upper left')

    ## position the labels appearing on the x axis
    x_positions = [np.log10(B) for B in x_labels]
    plt.xticks(x_positions, [format_B_label(B) for B in x_labels])

    ## plot the boxplots for A and B with different colors and offset in log space
    x_offset = 0.04
    width = 0.08
    draw_plot_at_positions(data_50, [x - x_offset for x in x_positions], width, dark_colors[2], light_colors[2])
    draw_plot_at_positions(data_100, [x + x_offset for x in x_positions], width, dark_colors[3], light_colors[3])

    plt.savefig("figures/Df_unif_varying_B.pdf")
    
def plot_Df_mignotte():
    ## recover the data in a convenient form
    data = {}
    print("Loading data for D(f) for Mignotte polynomials...")
    for x in list_input_Df_mignotte:
        degree = x[0]
        B = x[1]
        nb_tests = x[2]
        print("(degree,B) = ", degree,B)
        load("data/Df_"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+"_mignotte.sage")
        list_Df = []
        for y in res:
            list_Df += [y[0]]
        data[(degree,B)] = list_Df
    print("...done")

    list_degree = list(set([x[0] for x in list_input_Df_mignotte]))
    list_degree.sort()
    list_B = list(set([x[1] for x in list_input_Df_mignotte]))
    list_B.sort()

    ## create the plot
    plt.close()
        
    ## axes name, labels, legend
    fig, ax = plt.subplots()
    ax.set_ylabel('D(f)')
    ax.set_yscale('log') ## log scale on the y axis
    ax.set_xlabel('Degree of f')
    
    plt.xticks(range(len(list_degree)),list_degree)
    
    
    plt.legend(handles = [
        Patch(facecolor = light_colors[i%len(light_colors)], edgecolor=dark_colors[i%len(dark_colors)], label = "B = "+str(B)) for (i,B) in enumerate(list_B)], loc = 'upper left')
    
    
    ## plot the data
    for (i,B) in enumerate(list_B):
        y_data = [data[(degree,B)] for degree in list_degree if (degree,B) in data]
        draw_plot(y_data, 0.15*(i-len(list_B)/2), 0.15, dark_colors[i%len(dark_colors)], light_colors[i%len(light_colors)])

    ## save the plot
    plt.savefig("figures/Df_mignotte.pdf")


## TODO Uncomment below to create those plots (about D(f))
#plot_Df_unif_varying_degree()
#plot_Df_unif_varying_B()  
#plot_Df_mignotte()
    
#############################################
## Plots for primality of the discriminant ##
#############################################

def plot_prime_disc():
    ## recover the data in a convenient form
    data = {}
    print("Loading data for testing primality of discriminant...")
    for (degree,B,nb_tests) in list_input_prime_disc:
        print("(degree,B) = ", degree,B)
        load("data/prime_disc_"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+"_unif.sage")
        observed_number_of_prime = sum([x[0] for x in res])
        expected_number_of_prime = sum([x[1] for x in res])
        data[(degree,B)] = observed_number_of_prime/expected_number_of_prime
    print("...done")

    list_degree = list(set([x[0] for x in list_input_prime_disc]))
    list_degree.sort()
    list_B = list(set([x[1] for x in list_input_prime_disc]))
    list_B.sort()

    ## create the plot
    plt.close()
        
    ## axes name, labels, legend
    fig, ax = plt.subplots()
    ax.set_ylabel('observed value / expected value')
    ax.set_xlabel('Degree of f')
    
    plt.legend(handles = [
        Line2D([0],[0],color=dark_colors[i%len(dark_colors)], marker=list_markers[i%len(list_markers)], label = "B = "+str(B)) for (i,B) in enumerate(list_B)
    ], loc = 'upper left')
    
    ## plot the data
    for (i,B) in enumerate(list_B):
        x_data = [degree for degree in list_degree if (degree,B) in data]
        y_data = [data[(degree,B)] for degree in x_data]
        plt.plot(x_data,y_data, color=dark_colors[i%len(dark_colors)], marker=list_markers[i%len(list_markers)])

    ## save the plot
    plt.savefig("figures/prime_disc.pdf")
    
## TODO Uncomment below to create this plot (about primality of discriminant)
#plot_prime_disc()
    
#############################################
## Plots for norm of the conductor ideal   ##
#############################################

def plot_conductor():
    ## recover the data in a convenient form
    data = {}
    print("Loading data for norm of conductor ideals...")
    for (degree,B,nb_tests, max_time) in list_input_conductor:
        print("(degree,B) = ", degree,B)
        failed_input = 0
        load("data/norm_conductor_"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+"_unif.sage")
        list_conductors = []
        for x in res:
            if x[0] is None:
                failed_input += 1
            else:
                list_conductors += [x[0]]
        if np.median(list_conductors) != 1:
            print("\nfor (degree,B) = ", (degree,B), " the median is not 1 but ", np.median(list_conductors),"\n")
        data[(degree,B)] = list_conductors
        if failed_input > 0:
            print("\nFor (degree,B) = ", (degree,B), " there were ", Rshort(100*failed_input/nb_tests),"% failed inputs\n")
    print("...done")

    list_degree = list(set([x[0] for x in list_input_conductor]))
    list_degree.sort()
    list_B = list(set([x[1] for x in list_input_conductor]))
    list_B.sort()

    ## create the plot
    plt.close()
        
    ## axes name, labels, legend
    fig, ax = plt.subplots()
    ax.set_ylabel('normalized norm of the conductor ideal')
    ax.set_xlabel('Degree of f')
    
    plt.xticks(range(len(list_degree)),list_degree)
    
    
    plt.legend(handles = [
        Patch(facecolor = light_colors[i%len(light_colors)], edgecolor=dark_colors[i%len(dark_colors)], label = "B = "+str(B)) for (i,B) in enumerate(list_B)], loc = 'upper right')
    
    
    ## plot the data
    for (i,B) in enumerate(list_B):
        y_data = [data[(degree,B)] for degree in list_degree if (degree,B) in data]
        draw_plot(y_data, 0.15*(i-len(list_B)/2), 0.15, dark_colors[i%len(dark_colors)], light_colors[i%len(light_colors)])

    ## save the plot
    plt.savefig("figures/norm_conductor.pdf")
    
## TODO Uncomment below to create this plot (about the norm of the conductor ideal)
#plot_conductor()
    
######################################
## Plots for the Vandermonde norm   ##
######################################

def plot_Vandermonde(sampling_method ="unif", only_B_min = False):
    ## if only_B_min = True, plots the result only for the smallest value of B (i.e., B = 10 in the unif case and B = 100 in the mignotte case)
    ## in the mignotte case when only_B_min = True, the y axis is *not* in log-scale (in all the other cases it is in log-scale)
    
    ## recover the data in a convenient form
    data = {}
    print("Loading data for Vandermonde with sampling "+sampling_method+"...")
    if sampling_method == "unif":
        list_input = list_input_Vandermonde_unif
    elif sampling_method == "mignotte":
        list_input = list_input_Vandermonde_mignotte
    else:
        print("invalid sampling_method in plot_Vandermonde")
        return
    for y in list_input:
        degree = y[0]
        B = y[1]
        nb_tests = y[2]
        print("(degree,B) = ", degree,B)
        filename = "data/Vandermonde_"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+"_"+sampling_method+".sage"
        load(filename)
        list_Vandermonde_norms = []
        for x in res:
            if only_B_min and sampling_method == "mignotte":
                list_Vandermonde_norms += [x[-1]]
            else:
                list_Vandermonde_norms += [log(x[-1],10)]
        data[(degree,B)] = list_Vandermonde_norms
    print("...done")

    list_degree = list(set([x[0] for x in list_input]))
    list_degree.sort()
    list_B = list(set([x[1] for x in list_input]))
    list_B.sort()

    ## create the plot
    plt.close()
        
    ## axes name, labels, legend
    fig, ax = plt.subplots()
    if only_B_min and sampling_method == "mignotte":
        ax.set_ylabel(r'$\|V_f\|_F$')    
    else:
        ax.set_ylabel(r'$\log_{10}(\|V_f\|_F)$')
    #ax.set_yscale('log') ## log scale on the y axis
    ax.set_xlabel('Degree of f')
    
    plt.xticks(range(len(list_degree)),list_degree)
    
    if only_B_min:
        plt.legend(handles = [
            Patch(facecolor = light_colors[0], edgecolor=dark_colors[0], label = "B = "+str(list_B[0]))], loc = 'upper left')
    else:
        plt.legend(handles = [
            Patch(facecolor = light_colors[i%len(light_colors)], edgecolor=dark_colors[i%len(dark_colors)], label = "B = "+str(B)) for (i,B) in enumerate(list_B)], loc = 'upper left')
    
    
    ## plot the data
    if only_B_min:
        B = list_B[0]
        y_data = [data[(degree,B)] for degree in list_degree if (degree,B) in data]
        draw_plot(y_data, 0, 0.15, dark_colors[0], light_colors[0])
    else:
        for (i,B) in enumerate(list_B):
            y_data = [data[(degree,B)] for degree in list_degree if (degree,B) in data]
            draw_plot(y_data, 0.15*(i-len(list_B)/2), 0.15, dark_colors[i%len(dark_colors)], light_colors[i%len(light_colors)])

    ## save the plot
    if only_B_min:
        plt.savefig("figures/Vandermonde_B_"+str(B)+"_"+sampling_method+".pdf")
    else:
        plt.savefig("figures/Vandermonde_"+sampling_method+".pdf")

## TODO Uncomment below to create those plots (about the Frobenius norm of the Vandermonde matrix)
   
## We use latex in the plots below to have nicer font on the axes of the plots. If latex is not installed on the machine, this may raise a RuntimeError. In this case, we plot the figures without using latex (the name of the axes will be slightly less readable) 
    
#try:
#    plt.rcParams['text.usetex'] = True # to enable latex font
#    plot_Vandermonde(sampling_method ="unif", only_B_min = False)
#    plot_Vandermonde(sampling_method ="unif", only_B_min = True)
#    plot_Vandermonde(sampling_method ="mignotte", only_B_min = False)
#    plot_Vandermonde(sampling_method ="mignotte", only_B_min = True)
#    plt.rcParams.update(plt.rcParamsDefault) # restore matplotlib parameters to default
#except RuntimeError:
#    plt.rcParams.update(plt.rcParamsDefault) # restore matplotlib parameters to default
#    plot_Vandermonde(sampling_method ="unif", only_B_min = False)
#    plot_Vandermonde(sampling_method ="unif", only_B_min = True)
#    plot_Vandermonde(sampling_method ="mignotte", only_B_min = False)
#    plot_Vandermonde(sampling_method ="mignotte", only_B_min = True)

######################################################
## Plots for the minimum norm of H_y = Vf^T * Σ(y) ##
######################################################

## Note: norm_upperbound function is defined in test_min_norm_Hy.sage
## which is loaded via run_tests.sage (line 32)

def plot_min_norm_Hy(show_upperbound=False, only_B_min=False,
                     vary_B=False, fixed_degrees=None, per_degree=False):
    """
    Plot the minimum infinity norm of H_y matrices.

    Args:
        show_upperbound: If True, display upperbound asymptotic lines.
        only_B_min:      If True, only plot data for the smallest B value
                         (ignored when vary_B=True).
        vary_B:          If True, fix degree(s) and let B grow along the x-axis.
                         If False (default), fix B values and let degree grow.
        fixed_degrees:   Used only when vary_B=True. List of degrees to include
                         (e.g. [5, 6, 7]). If None, all available degrees are used.
        per_degree:      Used only when vary_B=True. If True, generate one separate
                         figure per degree (degree fixed, B grows). Each figure is
                         saved as min_norm_Hy_varying_B_deg{d}.pdf (or with
                         _with_upperbound suffix). If False (default), all degrees
                         are plotted together in a single figure.
    """
    ## recover the data in a convenient form
    data = {}
    print("Loading data for minimum norm of H_y...")
    for (degree, B, nb_tests, precision, zoom) in list_input_min_norm_Hy:
        print("(degree,B) = ", degree, B)
        filename = "data/norm_Hy"+str(degree)+"_"+str(B)+"_"+str(nb_tests)+".sage"
        try:
            load(filename)
            # res should be a list of norm_max_Hy values
            data[(degree, B)] = res
        except:
            print("  Warning: Could not load", filename)
            continue
    print("...done")

    list_degree = list(set([x[0] for x in list_input_min_norm_Hy]))
    list_degree.sort()
    list_B = list(set([x[1] for x in list_input_min_norm_Hy]))
    list_B.sort()

    ## ------------------------------------------------------------------ ##
    ##  Mode A: vary_B=True  →  x-axis = B, one series per fixed degree   ##
    ## ------------------------------------------------------------------ ##
    if vary_B:
        ## determine which degrees to plot
        all_degrees = list(list_degree)  # already sorted
        if fixed_degrees is None:
            plot_degrees = all_degrees
        else:
            plot_degrees = sorted([d for d in fixed_degrees if d in all_degrees])

        ## collect B values present for at least one of the requested degrees
        all_B_vary = sorted(list(set(
            [B for (d, B) in data.keys() if d in plot_degrees]
        )))

        if not all_B_vary:
            print("No data found for the requested degrees (vary_B mode).")
            return

        from matplotlib.ticker import FuncFormatter
        def log2_formatter_vB(y, pos):
            if y > 0:
                exponent = np.log2(y)
                if abs(exponent - round(exponent)) < 0.01:
                    return r'$2^{%d}$' % round(exponent)
                else:
                    return r'$2^{%.1f}$' % exponent
            return ''

        ## ------------------------------------------------------------------ ##
        ##  per_degree=True: one separate figure per degree                   ##
        ## ------------------------------------------------------------------ ##
        if per_degree:
            for (i_deg, degree) in enumerate(plot_degrees):
                ## B values available for this specific degree
                valid_B = sorted([B for B in all_B_vary if (degree, B) in data])
                if not valid_B:
                    print("  No data for degree", degree, "– skipping.")
                    continue

                plt.close()
                fig, ax = plt.subplots(figsize=(float(10), float(6)))

                ax.set_xlabel('Bound B on the coefficients')
                ax.set_ylabel(r'$\|V_f^\top \cdot T(y)\|_{max}$')
                # ax.set_title('degree = %d, varying B' % degree)
                ax.set_yscale('log', base=2)
                ax.yaxis.set_major_formatter(FuncFormatter(log2_formatter_vB))

                ## x-axis: log-scaled positions, labelled with B values
                x_positions = [np.log2(B) for B in valid_B]
                plt.xticks(x_positions,
                           [format_B_label(B) for B in valid_B])#, ha='right')

                edge_col = dark_colors[0]
                fill_col = light_colors[0]
                y_data   = [data[(degree, B)] for B in valid_B]
                pos      = [float(np.log2(B)) for B in valid_B]

                ax.boxplot(y_data,
                           positions=pos,
                           widths=float(0.25),
                           patch_artist=True,
                           manage_ticks=False,
                           showfliers=False,
                           whis=(10, 90),
                           boxprops=dict(facecolor=fill_col, edgecolor=edge_col,
                                         linewidth=float(2), alpha=float(0.7)),
                           medianprops=dict(color=edge_col, linewidth=float(2)),
                           whiskerprops=dict(color=edge_col, linewidth=float(1.5)),
                           capprops=dict(color=edge_col, linewidth=float(1.5)),
                           flierprops=dict(marker='o', markerfacecolor=edge_col,
                                           markersize=float(4), alpha=float(0.5)))

                legend_handles = [
                    Patch(facecolor=fill_col, edgecolor=edge_col,
                          label='degree = %d' % degree)
                ]

                ## optional upper-bound curve
                if show_upperbound:
                    ub_values = [float(norm_upperbound(degree, B)) for B in valid_B]
                    ax.plot(pos, ub_values,
                            color=edge_col,
                            linestyle='--', linewidth=float(2),
                            marker='*', markersize=float(10), alpha=float(0.8))
                    legend_handles.append(
                        Line2D([0], [0], color=edge_col, linestyle='--',
                               marker='*', markersize=float(8),
                               label=r'UB: $16m^5 B^5$ (m = %d)' % degree)
                    )

                ax.legend(handles=legend_handles, loc='upper left', framealpha=float(0.9))
                ax.grid(True, alpha=float(0.3), linestyle='--')

                if show_upperbound:
                    filename = "figures/min_norm_Hy_varying_B_deg%d_with_upperbound.pdf" % degree
                else:
                    filename = "figures/min_norm_Hy_varying_B_deg%d.pdf" % degree

                plt.savefig(filename, dpi=int(300), bbox_inches='tight')
                print("Plot saved to", filename)
            return

        ## ------------------------------------------------------------------ ##
        ##  per_degree=False (default): all degrees in one figure             ##
        ## ------------------------------------------------------------------ ##
        ## create the plot
        plt.close()
        fig, ax = plt.subplots(figsize=(float(10), float(6)))

        ax.set_xlabel('Bound B on the coefficients')
        ax.set_ylabel(r'$\|V_f^\top \cdot T(y)\|_{max}$')
        ax.set_yscale('log', base=2)
        ax.yaxis.set_major_formatter(FuncFormatter(log2_formatter_vB))

        ## log-scaled x positions, labelled with actual B values
        x_positions = [np.log2(B) for B in all_B_vary]
        plt.xticks(x_positions,
                   [format_B_label(B) for B in all_B_vary], rotation=30, ha='right')

        ## legend: one patch per degree
        legend_handles = [
            Patch(facecolor=light_colors[i % len(light_colors)],
                  edgecolor=dark_colors[i % len(dark_colors)],
                  label="degree = " + str(d))
            for (i, d) in enumerate(plot_degrees)
        ]

        ## box width / offset
        box_width_vB = float(0.6) / max(len(plot_degrees), 1)
        actual_box_width_vB = box_width_vB * 0.65

        for (i, degree) in enumerate(plot_degrees):
            color_idx = i % len(dark_colors)
            edge_col  = dark_colors[color_idx]
            fill_col  = light_colors[color_idx]

            valid_B = [B for B in all_B_vary if (degree, B) in data]
            y_data  = [data[(degree, B)] for B in valid_B]
            offset_val = box_width_vB * (i - len(plot_degrees) / 2.0)
            pos = [float(np.log2(B) + offset_val) for B in valid_B]

            if y_data:
                ax.boxplot(y_data,
                           positions=pos,
                           widths=actual_box_width_vB,
                           patch_artist=True,
                           manage_ticks=False,
                           showfliers=False,
                           whis=(10, 90),
                           boxprops=dict(facecolor=fill_col, edgecolor=edge_col,
                                         linewidth=float(2), alpha=float(0.7)),
                           medianprops=dict(color=edge_col, linewidth=float(2)),
                           whiskerprops=dict(color=edge_col, linewidth=float(1.5)),
                           capprops=dict(color=edge_col, linewidth=float(1.5)),
                           flierprops=dict(marker='o', markerfacecolor=edge_col,
                                           markersize=float(4), alpha=float(0.5)))

            ## optional upper-bound curve
            if show_upperbound:
                ub_values = [float(norm_upperbound(degree, B)) for B in valid_B]
                ub_pos    = [float(np.log2(B) + offset_val) for B in valid_B]
                ax.plot(ub_pos, ub_values,
                        color=edge_col,
                        linestyle='--', linewidth=float(2),
                        marker='*', markersize=float(10), alpha=float(0.8))
                legend_handles.append(
                    Line2D([0], [0], color=edge_col, linestyle='--',
                           marker='*', markersize=float(8),
                           label=r'UB: $16m^5 B^5$ (m = %d)' % degree)
                )

        ax.legend(handles=legend_handles, loc='upper left', framealpha=float(0.9))
        ax.grid(True, alpha=float(0.3), linestyle='--')

        ## build output filename
        if plot_degrees == all_degrees:
            deg_str = "all"
        else:
            deg_str = "_".join([str(d) for d in plot_degrees])
        if show_upperbound:
            filename = "figures/min_norm_Hy_varying_B_deg"+deg_str+"_with_upperbound.pdf"
        else:
            filename = "figures/min_norm_Hy_varying_B_deg"+deg_str+".pdf"

        plt.savefig(filename, dpi=int(300), bbox_inches='tight')
        print("Plot saved to", filename)
        return

    ## ------------------------------------------------------------------ ##
    ##  Mode B: vary_B=False (default) → x-axis = degree, one box per B  ##
    ## ------------------------------------------------------------------ ##
    if only_B_min:
        list_B = [list_B[0]]

    ## create the plot
    plt.close()
    
    ## Set up figure with publication-quality style
    fig, ax = plt.subplots(figsize=(float(10), float(6)))
    
    ## axes name
    ax.set_xlabel('Degree of f')
    ax.set_ylabel(r'$\|V_f^\top \cdot T(y)\|_{max}$')
    
    ## Set y-axis to log scale (base 2)
    ax.set_yscale('log', base=2)
    
    ## Format y-axis labels as 2^y
    from matplotlib.ticker import FuncFormatter
    def log2_formatter(y, pos):
        """Format y-axis labels as 2^y"""
        if y > 0:
            exponent = np.log2(y)
            if abs(exponent - round(exponent)) < 0.01:  # Close to an integer power
                return r'$2^{%d}$' % round(exponent)
            else:
                return r'$2^{%.1f}$' % exponent
        return ''
    ax.yaxis.set_major_formatter(FuncFormatter(log2_formatter))
    
    ## Format x-axis with log2 spacing
    # Use log2(degree) as x positions
    x_positions = [np.log2(d) for d in list_degree]
    x_labels = [str(d) for d in list_degree]
    
    plt.xticks(x_positions, x_labels)
    
    ## Create legend for box plots
    box_handles = []
    if only_B_min:
        box_handles = [
            Patch(facecolor=light_colors[0], edgecolor=dark_colors[0], 
                  label="B = "+str(list_B[0]))
        ]
    else:
        box_handles = [
            Patch(facecolor=light_colors[i%len(light_colors)], 
                  edgecolor=dark_colors[i%len(dark_colors)], 
                  label="B = "+str(B)) 
            for (i, B) in enumerate(list_B)
        ]
    
    ## plot the box plot data
    box_width = float(0.1)  # Base width
    actual_box_width = box_width * 0.8  # Actual rendered width
    
    for (i, B) in enumerate(list_B):
        # Only include degrees that have data for this B value
        valid_degrees = [degree for degree in list_degree if (degree, B) in data]
        y_data = [data[(degree, B)] for degree in valid_degrees]
        
        # Use log2(degree) as positions instead of indices
        positions_list = [np.log2(degree) for degree in valid_degrees]
        
        if y_data:  # Only plot if there's data
            if only_B_min:
                offset_val = float(0)
            else:
                offset_val = box_width * (i - len(list_B)/2.0)
           
            # Calculate actual positions with offset
            pos = [float(p + offset_val) for p in positions_list]
            
            # Set colors
            color_idx = 0 if only_B_min else i % len(dark_colors)
            edge_col = dark_colors[color_idx]
            fill_col = light_colors[color_idx]
            
            # Plot boxplot with proper styling matching reference
            bp = ax.boxplot(y_data, positions=pos, widths=actual_box_width, 
                           patch_artist=True, manage_ticks=False, 
                           showfliers=False, whis=(10,90),
                           boxprops=dict(facecolor=fill_col, edgecolor=edge_col, 
                                       linewidth=float(2), alpha=float(0.7)),
                           medianprops=dict(color=edge_col, linewidth=float(2)),
                           whiskerprops=dict(color=edge_col, linewidth=float(1.5)),
                           capprops=dict(color=edge_col, linewidth=float(1.5)),
                           flierprops=dict(marker='o', markerfacecolor=edge_col, 
                                          markersize=float(4), alpha=float(0.5)))
    
    ## Plot upperbound lines (if requested)
    line_handles = []
    if show_upperbound:
        for (i, B) in enumerate(list_B):
            # Only include degrees that have data for this B value
            valid_degrees = [degree for degree in list_degree if (degree, B) in data]
            
            if not valid_degrees:
                continue
            
            # Calculate upper bounds only for degrees with actual data
            ub_values = [float(norm_upperbound(d, B)) for d in valid_degrees]
            
            # Fit asymptotic formula in log-space
            log_degrees = np.log2(valid_degrees)
            log_ub = np.log2(ub_values)
            coeffs = np.polyfit(log_degrees, log_ub, 1)
            slope = float(coeffs[0])
            intercept = float(coeffs[1])
            
            # Create formula string for legend
            formula_str = r'$16m^5 B^5$ (B = %d)' %B
            
            # Calculate offset for this B value (same as box plot offset)
            if only_B_min:
                offset_val = float(0)
            else:
                offset_val = box_width * (i - len(list_B)/2.0)
            
            # Calculate x positions for degrees with data
            valid_x_positions = [np.log2(d) for d in valid_degrees]
            ub_x_positions = [float(x + offset_val) for x in valid_x_positions]
            
            # Plot upperbound points only (no connecting line, only star marker)
            # Ensure all parameters are Python native types
            line = ax.plot(ub_x_positions, ub_values,
                          color=dark_colors[i%len(dark_colors)],
                          linestyle='--', linewidth=float(2),
                          marker='*', markersize=float(12), alpha=float(0.8))
            line_handles.append((line[0], f'UB (B={B}): {formula_str}'))
    
    ## Combine legend handles
    all_handles = box_handles
    all_labels = [h.get_label() for h in box_handles]
    if show_upperbound and line_handles:
        all_handles += [h[0] for h in line_handles]
        all_labels += [h[1] for h in line_handles]
    
    ax.legend(all_handles, all_labels, loc='upper left', framealpha=float(0.9), ncol=2)
    ax.grid(True, alpha=float(0.3), linestyle='--')
    
    ## Expand y-axis upper limit to prevent legend from overlapping data
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax ) #* float(2**20)
    
    ## save the plot
    if show_upperbound:
        if only_B_min:
            filename = "figures/min_norm_Hy_B_"+str(list_B[0])+"_with_upperbound.pdf"
        else:
            filename = "figures/min_norm_Hy_with_upperbound.pdf"
    else:
        if only_B_min:
            filename = "figures/min_norm_Hy_B_"+str(list_B[0])+".pdf"
        else:
            filename = "figures/min_norm_Hy.pdf"
    
    plt.savefig(filename, dpi=int(300), bbox_inches='tight')
    print("Plot saved to", filename)

## TODO Uncomment below to create plots for minimum norm of H_y
plot_min_norm_Hy(show_upperbound=False, only_B_min=False)
plot_min_norm_Hy(show_upperbound=True, only_B_min=False)  # All B values with upperbound
plot_min_norm_Hy(show_upperbound=False, only_B_min=True)
plot_min_norm_Hy(show_upperbound=True, only_B_min=True)
## vary_B mode: fix degree(s), let B grow along x-axis (all in one figure)
#plot_min_norm_Hy(vary_B=True)                                                               # all degrees, no upper bound
#plot_min_norm_Hy(vary_B=True, show_upperbound=True)                                         # all degrees, with upper bound
#plot_min_norm_Hy(vary_B=True, fixed_degrees=[5, 6, 7, 8])                                   # select specific degrees
#plot_min_norm_Hy(vary_B=True, fixed_degrees=[5], show_upperbound=True)                      # single degree with UB

## per_degree mode: one separate figure per degree, degree fixed, B grows along x-axis
plot_min_norm_Hy(vary_B=True, per_degree=True, show_upperbound=True)                        # all degrees, each with UB
