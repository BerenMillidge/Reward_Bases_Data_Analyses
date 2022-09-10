import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# to do the fucking manual legend WHY WHY WHY FUCKING RAFAL ARGH
import matplotlib.patches as mpatch



neuronlist = ["359", "360", "361", "362", "363", "364", "365", "366", "367", "368", "369", "370", "371", "372", "373", "374", "375", "376","377"]

print(len(neuronlist))
# right let's try to get this stuff in the right format
def plot_lm_juices():
    coeflist = np.load("data/banana_juice_reward_type_coeff_list.npy")
    pslist = np.load("data/banana_juice_reward_type_ps_list.npy")
    #coeflist = np.load("lm_v1_coeffs_list.npy")
    #pslist = np.load("lm_v1_ps_list.npy")
    print(coeflist.shape)
    print(pslist.shape)
    coef_vals = []
    class_vals = []
    significants = []
    hue_vals = []
    for i in range(len(neuronlist)):
        # go get all the first ps then second ps
        coef = coeflist[2,i]
        pval = pslist[2,i]
        coef_vals.append(coef)
        if pval <= 0.05:
            significants.append(True)
        else:
            significants.append(False)
        hue_vals.append("Food")
        class_vals.append(neuronlist[i])
    for i in range(len(neuronlist)):
        # go get all the first ps then second ps
        coef = coeflist[3,i]
        pval = pslist[3,i]
        coef_vals.append(coef)
        if pval <= 0.05:
            significants.append(True)
        else:
            significants.append(False)
        hue_vals.append("Juice")
        class_vals.append(neuronlist[i])

    chi2_df = pd.DataFrame({'survived':hue_vals,
                        'class':class_vals,
                        'Frequency':coef_vals,
                        'sig':significants})

    plot(df_freq=chi2_df,x_var='class',hue_var="survived",title='Food/Juice Coefficients',sname="banana_juice_reward_type.jpg")
    
    
    
def plot_lm_value_separate_identity():
    coeflist = np.load("data/no_cross_gambles_lm_v1_reward_time_coeffs_list_both_monkeys.npy")
    pslist = np.load("data/no_cross_gambles_lm_v1_reward_time_ps_list_both_monkeys.npy")
    #coeflist = np.load("lm_v1_coeffs_list.npy")
    #pslist = np.load("lm_v1_ps_list.npy")
    print(coeflist.shape)
    print(pslist.shape)
    coef_vals = []
    class_vals = []
    significants = []
    labels = []
    for i in range(len(neuronlist)):
        # go get all the first ps then second ps
        coef = coeflist[1,i]
        pval = pslist[1,i]
        coef_vals.append(coef)
        if pval <= 0.05:
            significants.append(True)
        else:
            significants.append(False)
        labels.append("Value")
        class_vals.append(neuronlist[i])
    for i in range(len(neuronlist)):
        # go get all the first ps then second ps
        coef = coeflist[2,i]
        pval = pslist[2,i]
        coef_vals.append(coef)
        if pval <= 0.05:
            significants.append(True)
        else:
            significants.append(False)
        labels.append("Juice 1")
        class_vals.append(neuronlist[i])
    for i in range(len(neuronlist)):
        # go get all the first ps then second ps
        coef = coeflist[3,i]
        pval = pslist[3,i]
        coef_vals.append(coef)
        if pval <= 0.05:
            significants.append(True)
        else:
            significants.append(False)
        labels.append("Juice 2")
        class_vals.append(neuronlist[i])


    chi2_df = pd.DataFrame({'labels':labels,
                        'class':class_vals,
                        'Frequency':coef_vals,
                        'sig':significants})

    plot(df_freq=chi2_df,x_var='class',labels="labels",title='Value + Juice1 + Juice2 coefficients',sname="no_cross_gambles_reward_time_value_identity_separate_coefficients_lm_both_monkeys.jpg")

# okay great -- let's look at how to class this

# create a the dataframe
#chi2_df = pd.DataFrame({'survived':[0,0,0,1,1,1],
#                        'class':['First','Second','Third','First','Second','Third'],
#                        'Frequency':[80,97,372,136,87,119],
#                        'sig':[True,True,True,True,True,True]})


# plot results
# yeah it's really annoying having to mess with this when you don't understand it
#how to get the colors going the other wau
def sort_by_value():
    pass

# let's do all this sorting. super irritating for marginal gain and I don't see anyway to do it but they of course think it's possible
# because apparently all plotting is possible in matlab LIES!
from copy import deepcopy

def plot_lm_value_identity(neuronlist):
    coeflist = np.load("data/banana_juice_bipolar_coeff_list_normalized.npy")
    pslist = np.load("data/banana_juice_bipolar_ps_list_normalized.npy")
    
    coef_to_idx = {}
    for i in range(len(neuronlist)):
        coef_to_idx[coeflist[2,i]] = i
    first_coeffs = deepcopy(coeflist[2,:])
    first_coeffs.sort()
    first_coeffs = list(reversed(first_coeffs))
    print(list(first_coeffs))
    # now let's get everything sorted
    
    sorted_coeflist = np.zeros_like(coeflist)
    sorted_pslist = np.zeros_like(pslist)
    sorted_neuronlist = ["" for i in range(len(neuronlist))]
    for i in range(len(neuronlist)):
        coef = first_coeffs[i]
        idx = coef_to_idx[coef]
        sorted_coeflist[:,i] = coeflist[:,idx]
        sorted_pslist[:,i] = pslist[:,idx]
        sorted_neuronlist[i] = neuronlist[idx]
        
    coeflist = sorted_coeflist
    pslist = sorted_pslist
    neuronlist = sorted_neuronlist
    
    #coeflist = np.load("lm_v1_coeffs_list.npy")
    #pslist = np.load("lm_v1_ps_list.npy")
    print(coeflist.shape)
    print(pslist.shape)
    coef_vals = []
    class_vals = []
    significants = []
    labels = []
    for i in range(len(neuronlist)):
        # go get all the first ps then second ps
        coef = coeflist[1,i]
        pval = pslist[1,i]
        coef_vals.append(coef)
        if pval <= 0.05:
            significants.append(True)
        else:
            significants.append(False)
        labels.append("Value")
        class_vals.append(neuronlist[i])
    for i in range(len(neuronlist)):
        # go get all the first ps then second ps
        coef = coeflist[2,i]
        pval = pslist[2,i]
        coef_vals.append(coef)
        if pval <= 0.05:
            significants.append(True)
        else:
            significants.append(False)
        if coef > 0 or coef <= 0:
            labels.append("Identity")
        else:
            labels.append("Banana Responsive")
        class_vals.append(neuronlist[i])
    #for i in range(len(neuronlist)):
    #    if pval <= 0.05:
    #        significants.append(True)
    #    else:
    #        significants.append(False)
        
    

    chi2_df = pd.DataFrame({'labels':labels, # value or identity
                        'class':class_vals, # neuronlist
                        'Frequency':coef_vals,
                        'sig':significants})

    plot(df_freq=chi2_df,xvalues='class',labels="labels",title='Value/Identity coefficients',sname="banana_juice_value_identity_no_label.png")


# like this is completely crazy, I have NO IDEA what this is doing at all
def plot(df_freq,xvalues,labels,title,dst_dir=None,sname="", nolabel=True):

    stdres_sig_sorted = df_freq['sig'] #df_freq.sort_values(hue_var)['sig']
    print(df_freq)
    
    
    #palette=["gray","orange","blue"]
    coeflist = df_freq["Frequency"].to_numpy()
    print(coeflist.shape)
    colormat = [[],[]]
    colormat = []
    colormat += ["gray" for i in range(19)]
    colormat += ["orange" if coeflist[18 + i] > 0 else "blue" for i in range(19)]
    print(colormat)
    print(type(colormat))
    
    barplot = sns.barplot(x=xvalues, y='Frequency', hue=labels,data=df_freq,saturation=0.8,palette=["gray","blue"])
    #barplot.set_xticklabels(barplot.get_xticklabels(), rotation=45,fontsize=14)
    #barplot.get_xticklabels()
    if nolabel:
        barplot.set_xticklabels([])
    barplot.set_title(title,fontsize=16)
    barplot.set_xlabel("Neuron",fontsize=16)
    barplot.set_ylabel("Coefficient",fontsize=16)
    #barplot.set_yticklabels(barplot.get_yticklabels(),fontsize=14)
    # assign legend first to prevent color overlap
    #plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.5), ncol=1)
    plt.yticks(fontsize=12)
    
    # manual legend creation. Fuck rafal. this is so fucking stupid
    grey_patch = mpatch.Patch(color="gray", label = 'Subjective Value')
    blue_patch = mpatch.Patch(color = 'blue', label='Identity (Juice)')
    orange_patch = mpatch.Patch(color='orange', label='Identity (Banana)')
    
    
    plt.legend(handles = [grey_patch, blue_patch, orange_patch],loc='upper left',fontsize=10)

    print(barplot.patches)

    for p, sig in zip(barplot.patches,stdres_sig_sorted):
        if sig == True:
            barplot.text(p.get_x() + p.get_width() / 2., p.get_height(),
                        '*', ha='center')
            #barplot.axvline(p.get_x() + p.get_width() / 2., lw=0.5)
    
    for i,bar in enumerate(barplot.patches):
        if i > 18:
            if coeflist[i] < 0:
                bar.set_color('orange')
            bar.set_alpha(0.8)




    plt.tight_layout()
    plt.savefig("figures/" + str(sname), format="png")
    plt.show()


if __name__ == '__main__':
    #plot_lm_juices()
    plot_lm_value_identity(neuronlist)
    #plot_lm_value_separate_identity()

# great so this is nice. they almost all respond strongly to identity.
# now I need to do this shuffling technique arbitrarily to obtain some kind of p value
#  but it's looking vaguely promising although a massive effect of value, so I have
# right so let's thinkg about how well it works.
# how should I actually do this...
# one thing I can do is plot how many are significant and then randomize
