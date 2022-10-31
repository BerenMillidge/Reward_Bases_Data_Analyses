# code for simulations of toy model task with weighted reward prediction errors

from cProfile import run
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import os

#5 = 0.9ml juice subj value 1
#1 = 1.5kg bananas 0.7
# 4 = 0.5ml juice 0.6
# 3 = 0.2ml juice 0.1
# 2  = 0.3kg bananas 0.05


neuronlist = ["359", "360", "361", "362", "363", "364", "365", "366", "367", "368", "369", "370", "371", "372", "373", "374", "375", "376","377"]


def run_simulation(N_runs = 200, lr = 0.1, N_seeds = 10, regression_type = "reward_type_only", double_update = False,use_subjective_value_reward = True):
    Vs_food = np.zeros(5)
    Vs_water = np.zeros(5)
    if use_subjective_value_reward:
        rs_food = np.array([0.7, 0.05,0.0,0.0,0.0])
        rs_water = np.array([0.0,0.0,0.1,0.6,1.0])
    else:
        # use reward_magnitude
        rs_food = np.array([1.5,0.3,0.0,0.0,0.0])
        rs_water = np.array([0.0,0.0,0.2,0.5,0.9])
    
    if regression_type == "subjective_value_reward_type":
        coeflist = np.load("data/banana_juice_reward_type_coeff_list_normalized.npy")
        juice_coeffs = coeflist[2,:]
        food_coeffs = coeflist[3,:]
        N = len(juice_coeffs)
    elif regression_type == "reward_type_only":
        coeflist = np.load("data/reward_type_only_coeff_list_normalized.npy")
        juice_coeffs = coeflist[1,:]
        food_coeffs = coeflist[2,:]
        N = len(juice_coeffs)
    elif regression_type == "bipolar_identity":
        coeflist = np.load("data/banana_juice_bipolar_coeff_list_normalized.npy")
        val_coeffs = coeflist[1,:]
        identity_coeffs = coeflist[2,:]
        N = len(identity_coeffs)
    elif regression_type == "reward_type_subjective_value":
        coeflist = np.load("data/reward_type_only_subjective_value_coeff_list_normalized.npy")
        juice_coeffs = coeflist[1,:]
        food_coeffs = coeflist[2,:]
        N = len(juice_coeffs)
    else:
        raise ValueError("regression type not recognized")

    Vsss_food = []
    Vsss_water = []

    for j in range(N_seeds):
        print("N seeds: ", N_seeds)
        Vss_food = []
        Vss_water = []
        Vs_food = np.zeros(5)
        Vs_water = np.zeros(5)
        for i in range(N_runs):
            idx = int(np.random.uniform(low=0,high=5))
            delta_food = rs_food[idx] - Vs_food[idx]
            delta_water = rs_water[idx] - Vs_water[idx]
            delta_V_food = 0
            delta_V_water = 0
            for n in range(N):
                if regression_type == "bipolar_identity":
                    delta_V_food += identity_coeffs[n] * delta_food
                    delta_V_water += identity_coeffs[n] * delta_water
                elif regression_type == "reward_type_only" or regression_type == "subjective_value_reward_type" or regression_type == "reward_type_subjective_value":
                    delta_V_food += food_coeffs[n] * delta_food
                    delta_V_water += juice_coeffs[n] * delta_water
                    if double_update:
                        delta_V_food += juice_coeffs[n] * delta_water
                        delta_V_water += food_coeffs[n] * delta_food
            Vs_food[idx] += lr * delta_V_food
            Vs_water[idx] += lr * delta_V_water
            Vss_food.append(deepcopy(Vs_food))
            Vss_water.append(deepcopy(Vs_water))
        Vss_food = np.array(Vss_food)
        Vss_water = np.array(Vss_water)
        Vsss_food.append(deepcopy(Vss_food))
        Vsss_water.append(deepcopy(Vss_water))

    Vsss_food = np.array(Vsss_food)
    Vsss_water = np.array(Vsss_water)
    print("vssss")
    print(Vsss_food.shape)
    print(Vsss_water.shape)
    mean_food = np.mean(Vsss_food, axis=0)
    mean_water = np.mean(Vsss_water, axis=0)
    std_food = np.std(Vsss_food, axis=0) / np.sqrt(N_seeds)
    std_water = np.std(Vsss_water, axis=0) / np.sqrt(N_seeds)
    # begin plot

    fig = plt.figure(figsize=(10,10))
    xs = np.array([i for i in range(N_runs)])
    print(xs.shape)
    print(mean_food.shape)
    print(std_food.shape)
    LINEWIDTH = 2
    plt.plot(xs, mean_food[:,0], label="Banana 1.5g", color="orange", linestyle="solid", linewidth = LINEWIDTH)
    plt.fill_between(xs, mean_food[:,0] - std_food[:,0], mean_food[:,0] + std_food[:,0], color="orange", alpha=0.5)
    plt.hlines(rs_food[0],xmin=0, xmax=N_runs, color="orange", linestyle="dashed", alpha=0.5, linewidth = LINEWIDTH)
    plt.plot(xs, mean_food[:,1], label="Banana 0.3g", color="orange", linestyle="dashed", linewidth = LINEWIDTH)
    plt.fill_between(xs, mean_food[:,1] - std_food[:,1], mean_food[:,1] + std_food[:,1], color="orange", alpha=0.3)
    plt.hlines(rs_food[1],xmin=0, xmax=N_runs, color="orange", linestyle="dashed", alpha=0.5, linewidth = LINEWIDTH)
    plt.plot(xs, mean_water[:,4], label="Juice 0.9ml", color="blue", linestyle="solid", linewidth = LINEWIDTH)
    plt.fill_between(xs, mean_water[:,4] - std_water[:,4], mean_water[:,4] + std_water[:,4], color="blue", alpha=0.5)
    plt.hlines(rs_water[4],xmin=0, xmax=N_runs, color="blue", linestyle="dashed", alpha=0.5, linewidth = LINEWIDTH)
    plt.plot(xs, mean_water[:,3], label="Juice 0.5ml", color="blue", linestyle="dashed", linewidth = LINEWIDTH)
    plt.fill_between(xs, mean_water[:,3] - std_water[:,3], mean_water[:,3] + std_water[:,3], color="blue", alpha=0.3)
    plt.hlines(rs_water[3],xmin=0, xmax=N_runs, color="blue", linestyle="dashed", alpha=0.5, linewidth = LINEWIDTH)
    plt.plot(xs, mean_water[:,2], label="Juice 0.2ml", color="blue", linestyle="dotted", linewidth = LINEWIDTH)
    plt.fill_between(xs, mean_water[:,2] - std_water[:,2], mean_water[:,2] + std_water[:,2], color="blue", alpha=0.1)
    plt.hlines(rs_water[2],xmin=0, xmax=N_runs,color="blue", linestyle="dashed", alpha=0.5, linewidth = LINEWIDTH)
    if use_subjective_value_reward:
        plt.ylabel("Subjective value function",fontsize=28)
    else:
        plt.ylabel("Value function",fontsize=28)
    plt.xlabel("Trial",fontsize=28)
    plt.title("Learning the value function",fontsize=28)
    plt.legend(loc="lower right", fontsize=25)
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    plt.tight_layout()
    fig.savefig("figures/value_function_learning_" + str(regression_type) + ".png", format="png")
    plt.show()

# this is also in the julia code!
def plot_coefficient_scatter(regression_type = "reward_type_only"):  
    if regression_type == "subjective_value_reward_type":
        coeflist = np.load("data/banana_juice_reward_type_coeff_list_normalized.npy")
        juice_coeffs = coeflist[2,:]
        food_coeffs = coeflist[3,:]
    elif regression_type == "reward_type_only":
        coeflist = np.load("data/reward_type_only_coeff_list_normalized.npy")
        juice_coeffs = coeflist[1,:]
        food_coeffs = coeflist[2,:]
    elif regression_type == "reward_type_subjective_value":
        coeflist = np.load("data/reward_type_only_subjective_value_coeff_list_normalized.npy")
        juice_coeffs = coeflist[1,:]
        food_coeffs = coeflist[2,:]
    else:
        raise ValueError("regression type not recognized")
    fig = plt.figure(figsize=(12,10))
    # print out values to determine which neurons are used
    print(neuronlist[0], neuronlist[9])
    print(juice_coeffs[0], food_coeffs[0])
    print(juice_coeffs[9], food_coeffs[9])
    plt.scatter(juice_coeffs, food_coeffs, color='black')
    plt.hlines(0,xmin=-1,xmax=2, color='gray', linestyle="dashed")
    plt.vlines(0,ymin=-1, ymax=2, color='gray', linestyle="dashed")
    plt.xlabel("Banana coefficient",fontsize=28)
    plt.ylabel("Juice coefficient",fontsize=28)
    plt.ylim(-0.5,2)
    plt.xlim(-0.5,2)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.title("Regression coefficients per neuron",fontsize=28)
    plt.savefig("figures/food_juice_scatterplot_" + str(regression_type) + ".png", format="png")
    plt.gca().set_aspect('equal') # for rafal's fucking aspect ratio requests
    plt.tight_layout()
    plt.show()
            
def plot_simulated_neuron_timecourse(regression_type="reward_type_only", no_xticks = True):
    FONT_SIZE = 30
    if regression_type == "subjective_value_reward_type":
        coeflist = np.load("data/banana_juice_reward_type_coeff_list_normalized.npy")
        juice_coeffs = coeflist[3,:]
        food_coeffs = coeflist[2,:]
    elif regression_type == "reward_type_only":
        coeflist = np.load("data/reward_type_only_coeff_list_normalized.npy")
        juice_coeffs = coeflist[1,:]
        food_coeffs = coeflist[2,:]
    elif regression_type == "reward_type_subjective_value":
        coeflist = np.load("data/reward_type_only_subjective_value_coeff_list_normalized.npy")
        juice_coeffs = coeflist[1,:]
        food_coeffs = coeflist[2,:]

    else:
        raise ValueError("regression type not recognized")
    #plot_coefficient_scatter()
    print("JUICE_COEFFS: ", juice_coeffs)
    print("FOOD COEFFS: ", food_coeffs)
    food_juice_diffs = juice_coeffs - food_coeffs
    print("JUICE_FOOD_DIFFS: ",food_juice_diffs)
    max_juice_idx =np.argmax(food_juice_diffs)
    min_juice_idx = np.argmin(food_juice_diffs)
    print(juice_coeffs[max_juice_idx], food_coeffs[max_juice_idx])
    print(juice_coeffs[min_juice_idx], food_coeffs[min_juice_idx])
    # let's plot the 'max juice neuron
    # let's do it for each condition to see what's up!
    xs = np.array([-200,0,200,400,600])
    rs_conditions = np.array([1.5,0.3,0.2,0.5,0.9])
    rs_food = np.array([1.5,0.3,0.0,0.0,0.0])
    rs_water = np.array([0.0,0.0,0.2,0.5,0.9])
    avg_food = np.mean(rs_conditions[0:2])
    avg_water = np.mean(rs_conditions[2:len(rs_conditions)])
    print(avg_food)
    print(avg_water)
    condition_labels = ["1.5g Banana", "0.3g Banana", "0.2ml Juice","0.5ml Juice","0.9ml Juice"]
    predicted_vals_min = np.zeros((5, len(xs)))
    predicted_vals_max = np.zeros((5, len(xs)))
    # so yeah we can generate vaguely desirable plots but obviously nothing interesting actually happens after but this is def dealbale
    # with a bit of hacking around the averages!
    linestyles = ["solid","dashed","dotted","dashed","solid"]
    colors = ["orange","orange","blue","blue","blue"]
    for i in range(5):
        if i < 2:
            avg_food_temp = avg_food
        else:
            avg_food_temp = 0
        if i >= 2:
            #avg_water_temp = 0.0 # a total numerical hack to fix this but it kind of works
            avg_water_temp = 0
        else:
            avg_water_temp = avg_water
        print("IDX: ", i)
        delta_min = (food_coeffs[min_juice_idx] * (rs_food[i] - avg_food_temp)) + (juice_coeffs[min_juice_idx] * (rs_water[i] - avg_water_temp))
        predicted_vals_min[i,2] = delta_min
        print("delta min: ", delta_min)
        delta_max = (food_coeffs[max_juice_idx] * (rs_food[i] - avg_food_temp)) + (juice_coeffs[max_juice_idx] * (rs_water[i] - avg_water_temp))
        predicted_vals_max[i,2] = delta_max
        print("delta max: ", delta_max)
        

    fig = plt.figure(figsize=(10,10))
    for i in range(5):
        plt.plot(xs,predicted_vals_min[i,:],color = colors[i], linestyle=linestyles[i],label = condition_labels[i])
    plt.legend(fontsize=24)
    plt.xlabel("Time (ms)", fontsize=FONT_SIZE)
    plt.ylabel("Simulated spike count", fontsize=FONT_SIZE)
    plt.vlines(0,ymin=np.amin(predicted_vals_min[:,2]),ymax=np.amax(predicted_vals_min[:,2]),color="gray", linestyle="dashed", label="Cue Onset")
    if no_xticks:
        plt.xticks([])
    else:
        plt.xticks(fontsize=FONT_SIZE)
    plt.yticks(fontsize=FONT_SIZE)
    plt.title("Simulated banana neuron",fontsize=FONT_SIZE)
    plt.tight_layout()
    fig.savefig("figures/simulated_banana_neuron_" + str(regression_type) + ".png", format="png")
    plt.show()

    fig = plt.figure(figsize=(10,10))
    for i in range(5):
        plt.plot(xs,predicted_vals_max[i,:],color = colors[i], linestyle=linestyles[i],label = condition_labels[i])
    plt.legend(fontsize=24)
    plt.xlabel("Time (ms)", fontsize=FONT_SIZE)
    plt.ylabel("Simulated spike count", fontsize=FONT_SIZE)
    plt.vlines(0,ymin=np.amin(predicted_vals_max[:,2]),ymax=np.amax(predicted_vals_max[:,2]),color="gray", linestyle="dashed", label="Cue Onset")
    if no_xticks:
        plt.xticks([])
    else:
        plt.xticks(fontsize=FONT_SIZE)
    plt.yticks(fontsize=FONT_SIZE)
    plt.title("Simulated juice neuron",fontsize=FONT_SIZE)
    plt.tight_layout()
    fig.savefig("figures/simulated_juice_neuron_" + str(regression_type) + ".png", format="png")
    plt.show()
        
        
if __name__ == '__main__':
    if not os.path.exists("figures"):
        os.makedirs("figures")
    N_runs = 50
    N_seeds = 1000
    regression_type = "reward_type_subjective_value"
    
    #plot_simulated_neuron_timecourse(regression_type=regression_type)
    #plot_coefficient_scatter(regression_type=regression_type) 
    run_simulation(N_runs = N_runs,N_seeds = N_seeds,double_update=False,regression_type=regression_type)



#coeflist = np.load("banana_juice_reward_type_coeff_list_normalized.npy")
#val_coeffs = coeflist[1,:]
#juice_coeffs = coeflist[2,:]
#food_coeffs = coeflist[3,:]
#reg_coeffs = coeflist[0,:]
#print(val_coeffs)
#print(reg_coeffs)

        
        
# put only 100-150 things
# simulation multipanel, results of simulation
# putting some schematics of how RB works
# useful to show scatterplot of coefficients and and b across neurons
# third plot illustrating responses to simulated neurons equivalent to real data
# plot for banana neuron, juice neuron, value neuron, plot timelines of prediction errors
# subtract average amoutn across five curves corresponding to five curves in the other plots and show how it looks like
# i.e. simulating learning with neurons which behave like what armin observed
# actually will need to get these binary conditions to see what's up
        
    

    


