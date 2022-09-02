# plotting functions for Fig1
using MAT
using Plots
using GLM
using DataFrames
using Distributions
using StatsBase
using HypothesisTests
using StatsPlots
using Statistics
using NPZ
using Plots.PlotMeasures
include("utils.jl")

ALL_SITUATIONS = [1,2,3,4,5,25]
neuronlist = [359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377]

juice_neuron = [359]
food_neuron = [368]
#food_neuron = [370]
bucket_start = -500
bucket_end = 1500
window_size = 100
standard_neuron = [366]

ALL_JUICE_SITUATIONS = [3,4,5]
ALL_BANANA_SITUATIONS = [1,2]
JUICE_SITUATIONS = [5]
BANANA_SITUATIONS = [1]
EXPERIMENT_SITUATIONS = [1,2,3,4,5]

GUIDEFONTSIZE=15
TICKFONTSIZE=10
LEGENDFONTSIZE=10
TITLEFONTSIZE=15

function mean_and_std_spikes()
    spikelist = []
    for neuron in neuronlist
        spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" * ".mat")
        for s in spiketimes
            push!(spikelist, length(s))
        end
    end
    return mean(spikelist), std(spikelist) / sqrt(length(neuronlist)), spikelist
end

function plot_average_spikecourse()
    bucket_start = -2500
    bucket_end = 8000
    spikes = bucket_count_spikes(stim_relative_times, 100,bucket_start, bucket_end)
    xs = collect(bucket_start:100:bucket_end)
    plot(xs, spikes)
end

function plot_overlapping_window_timecourses(window_size=100, use_frequency_labels = true)
    bucket_start = -500
    bucket_end = 1500
    window_size=window_size
    xs = collect(bucket_start:50:bucket_end-100)
    juice_means, juice_stds = count_plot_by_neuron(juice_neuron,bucket_start, bucket_end, window_size,true,true)
    food_means, food_stds = count_plot_by_neuron(food_neuron, bucket_start, bucket_end, window_size,true,true)
    standard_means, standard_stds = count_plot_by_neuron(standard_neuron, bucket_start, bucket_end, window_size,true,true)
    if use_frequency_labels
        window_size_seconds = window_size / 1000
        juice_means = juice_means ./ window_size_seconds
        food_means = food_means ./ window_size_seconds
    end
    plot(xs,juice_means[5][1,:], label="0.9ml juice",color=:blue, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE,legendfontsize=LEGENDFONTSIZE)
    plot!(xs,juice_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
    plot!(xs,juice_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
    plot!(xs,juice_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
    plot!(xs,juice_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
    vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
    xlabel!("Time (ms) after cue")
    if use_frequency_labels
        ylabel!("Firing rate within window (Hz)")
    else
        ylabel!("Spikes within window")
    end
    title!("Juice responsive neuron")
    savefig("figures/juice_neuron_all_conditions_overlapping_window_2.png")

    plot(xs,food_means[5][1,:], label="0.9ml juice",color=:blue, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE,legendfontsize=LEGENDFONTSIZE)
    plot!(xs,food_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
    plot!(xs,food_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
    plot!(xs,food_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
    plot!(xs,food_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
    vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
    xlabel!("Time (ms) after cue")
    if use_frequency_labels
        ylabel!("Firing rate within window (Hz)")
    else
        ylabel!("Spikes within window")
    end
    title!("Food responsive neuron")
    savefig("figures/food_neuron_all_conditions_overlapping_window_plot_2.png")

    plot(xs,standard_means[5][1,:], label="0.9ml juice",color=:blue, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE,legendfontsize=LEGENDFONTSIZE)
    plot!(xs,standard_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
    plot!(xs,standard_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
    plot!(xs,standard_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
    plot!(xs,standard_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
    vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
    xlabel!("Time (ms) after cue")
    if use_frequency_labels
        ylabel!("Firing rate within window (Hz)")
    else
        ylabel!("Spikes within window")
    end
    title!("Value responsive neuron")
    savefig("figures/value_only_neuron_all_conditions_overlapping_window_plot_2.png")
end

function timecourse_plots_non_overlapping()
    bucket_start = -500
    bucket_end = 1500
    window_size=100
    xs = collect(bucket_start:window_size:bucket_end)
    juice_means, juice_stds = count_plot_by_neuron(juice_neuron,bucket_start, bucket_end, window_size,true)
    food_means, food_stds = count_plot_by_neuron(food_neuron, bucket_start, bucket_end, window_size,true)
    standard_means, standard_stds = count_plot_by_neuron(standard_neuron, bucket_start, bucket_end, window_size,true)
    plot(xs[2:end-1],juice_means[5][1,:], label="0.9ml juice",color=:blue, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE, legendfontsize=LEGENDFONTSIZE)
    plot!(xs[2:end-1],juice_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
    plot!(xs[2:end-1],juice_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
    plot!(xs[2:end-1],juice_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
    plot!(xs[2:end-1],juice_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
    vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
    xlabel!("Time (ms) after cue")
    ylabel!("Spikes within window")
    title!("Juice Responsive Neuron")
    savefig("figuresjuice_neuron_all_conditions_window_plot_3.png")

    plot(xs[2:end-1],food_means[5][1,:], label="0.9ml juice",color=:blue, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE, legendfontsize=LEGENDFONTSIZE)
    plot!(xs[2:end-1],food_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
    plot!(xs[2:end-1],food_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
    plot!(xs[2:end-1],food_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
    plot!(xs[2:end-1],food_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
    vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
    xlabel!("Time (ms) after cue")
    ylabel!("Spikes within window")
    title!("Food Responsive Neuron")
    savefig("figuresfood_neuron_all_conditions_window_plot_3.png")

    plot(xs[2:end-1],standard_means[5][1,:], label="0.9ml juice",color=:blue, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE, legendfontsize=LEGENDFONTSIZE)
    plot!(xs[2:end-1],standard_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
    plot!(xs[2:end-1],standard_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
    plot!(xs[2:end-1],standard_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
    plot!(xs[2:end-1],standard_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
    vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
    xlabel!("Time (ms) after cue")
    ylabel!("Spikes within window")
    title!("Spikes of Value Responsive Neuron")
    savefig("figuresvalue_only_neuron_all_conditions_window_plot_2.png")
end

function run_lms()
    window_start = 100
    window_end = 400
    coeff_list, ts_list, ps_list = subjective_value_reward_type_coeffs(neuronlist, ALL_SITUATIONS, window_start,window_end, false, true)
    npzwrite("banana_juice_reward_type_coeff_list_normalized.npy", coeff_list)
    npzwrite("banana_juice_reward_type_ts_list_normalized.npy", ts_list)
    npzwrite("banana_juice_reward_type_ps_list_normalized.npy", ps_list)

    coeff_list, ts_list, ps_list = subjective_value_bipolar_reward(neuronlist, ALL_SITUATIONS,window_start,window_end,false, true)
    npzwrite("banana_juice_bipolar_coeff_list_normalized.npy", coeff_list)
    npzwrite("banana_juice_bipolar_ts_list_normalized.npy", ts_list)
    npzwrite("banana_juice_bipolar_ps_list_normalized.npy", ps_list)
end

function juice_neuron_response_histogram()
    window_start = 100
    window_end = 400
    spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0359" *".mat")
    count_list_juice_neuron_juice = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, JUICE_SITUATIONS,window_start,window_end,false)
    count_list_juice_neuron_food = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, BANANA_SITUATIONS, window_start,window_end,false)
    xs = collect(1:11)
    histogram(xs,count_list_juice_neuron_juice,nbins=10,alpha=0.5,label="Juice Response", guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE, legendfontsize=LEGENDFONTSIZE)
    histogram!(xs,count_list_juice_neuron_food,nbins=10, alpha=0.5, label="Food Response")
    xlabel!("Spike Count in Window")
    ylabel!("Number of Trials")
    title!("Juice Neuron Response to Food/Juice")
    savefig("figures/juice_neuron_response_histogram_2.png")
end


function all_neurons_condition_histogram(window_start = 100, window_end = 400,remove_labels = false,plot_frequency=false, plot_separate_histograms=false)
    if plot_frequency
        window_size = window_end - window_start
        divisor = window_size / 1000 # msec to sec
    else
        divisor = 1
    end
    spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0359" *".mat")
    situation_1 = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, [1],window_start,window_end,false) / divisor
    situation_2 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [2], window_start,window_end,false) / divisor
    situation_3 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [3], window_start,window_end,false) / divisor
    situation_4 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [4], window_start,window_end,false) / divisor
    situation_5 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [5], window_start,window_end,false) / divisor
    #histogram(xs,[situation_1,situation_2,situation_3,situation_4,situation_5],nbins=10,alpha=0.5,label=["1.5kg Banana" "0.3kg Banana" "0.2ml Juice" "0.5ml Juice" "0.9ml Juice"],layout=(5,1),xlims=(0,9),ylims=(0,6),xticks=xs)
    spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0368" *".mat")
    situation_2_1 = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, [1],window_start,window_end,false) / divisor
    situation_2_2 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [2],window_start,window_end,false) / divisor
    situation_2_3 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [3],window_start,window_end,false) / divisor
    situation_2_4 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [4],window_start,window_end,false) / divisor
    situation_2_5 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [5],window_start,window_end,false) / divisor
    if plot_frequency
        xs = collect(0:3:30)
    else
        xs = collect(0:1:10)
    end
    bar_width = (1/divisor) - 0.3
    xticks = (xs[1:end-1] .+ (0.5 / divisor), string.(xs))
    histogram_list = [situation_1,situation_2_1,situation_2, situation_2_2,situation_3,situation_2_3,situation_4,situation_2_4,situation_5,situation_2_5]
    label_list = ["1.5g Banana" "1.5g Banana" "0.3g Banana" "0.3g Banana" "0.2ml Juice" "0.2ml Juice" "0.5ml Juice" "0.5ml Juice" "0.9ml Juice" "0.9ml Juice"]
    color_list = [:orange :orange :orange :orange :blue :blue :blue :blue :blue :blue]
    if !plot_separate_histograms
        if remove_labels
            hist_plot = histogram(xs,[situation_1,situation_2_1,situation_2, situation_2_2,situation_3,situation_2_3,situation_4,situation_2_4,situation_5,situation_2_5],alpha=0.5,label=["1.5g Banana" "1.5g Banana" "0.3g Banana" "0.3g Banana" "0.2ml Juice" "0.2ml Juice" "0.5ml Juice" "0.5ml Juice" "0.9ml Juice" "0.9ml Juice"],layout=(5,2),size=(1000,800),ylims=(0,6), bins=xs,xticks=xticks, color = [:orange :orange :orange :orange :blue :blue :blue :blue :blue :blue],bar_width=bar_width, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE, legendfontsize=LEGENDFONTSIZE)
        else
            hist_plot = histogram(xs,[situation_1,situation_2_1,situation_2, situation_2_2,situation_3,situation_2_3,situation_4,situation_2_4,situation_5,situation_2_5],alpha=0.5,label=["1.5g Banana" "1.5g Banana" "0.3g Banana" "0.3g Banana" "0.2ml Juice" "0.2ml Juice" "0.5ml Juice" "0.5ml Juice" "0.9ml Juice" "0.9ml Juice"],layout=(5,2),xlims=(0,10),ylims=(0,6),xticks=xticks,size=(1000,800),xlabel="Spike Count in Window", ylabel="Number of Trials", color = [:orange :orange :orange :orange :blue :blue :blue :blue :blue :blue],bar_width=bar_width,bins = xs, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE, legendfontsize=LEGENDFONTSIZE)
        end
        savefig("figures/all_neurons_conditions_histogram_5.png")
    else
        N = length(histogram_list)
        for i in 1:N
            hist_plot = histogram(xs, histogram_list[i],label=label_list[i],color=color_list[i],ylims=(0,6), bins=xs,xticks=xticks,bar_width=bar_width, guidefontsize=GUIDEFONTSIZE, tickfontsize=15, legendfontsize=14,size=(500,125),yticks=[2,4,6])
            savefig("figures/all_neurons_condition_histogram_split_2_$i" * ".png")
        end
    end
end

function bar_chart_spike_proportions()
    window_start = 100
    window_end = 400
    mean_1s, mean_2s, stds_1s, stds_2s = neurons_spike_proportions(neuronlist, JUICE_SITUATIONS, BANANA_SITUATIONS,window_start,window_end)
    bar_xticks = (collect(1:L),string.(neuronlist)[1,:] )
    bar(mean_1s,yerr = stds_1s ./ L, legend = :none, xticks = bar_xticks, guidefontsize=GUIDEFONTSIZE, tickfontsize=TICKFONTSIZE, legendfontsize=LEGENDFONTSIZE)
    xlabel!("Neuron")
    ylabel!("Spike Proportion Juice")
    title!("Proportion of total spikes to juice or food reward")
    savefig("figures/proportions_bar_100_400_2.png")

    bar_xticks = (collect(1:L),string.(neuronlist)[1,:] )
    bar(mean_1s .- mean(mean_1s),yerr = stds_1s ./ L, legend = :none, xticks = bar_xticks)
    xlabel!("Neuron")
    ylabel!("Proportion Difference from mean response")
    title!("Proportion Differences to juice or food reward")
    savefig("figures/proportions_diff_100_400_2.png")
end

function relative_spikes_scatterplot_old(plot_firing_rate=false)
    window_start = 100
    window_end = 400
    window_size = window_end - window_start
    divisor = window_size / 1000 # msec to sec conversion
    counts_juice, counts_banana = neuron_means_std_spikes_in_window_by_situation(neuronlist,[JUICE_SITUATIONS, BANANA_SITUATIONS],window_start,window_end,false, true)
    counts_juice_means, counts_juice_stds = counts_juice
    counts_banana_means, counts_banana_stds = counts_banana
    if plot_firing_rate
        counts_juice_means = counts_juice_means / divisor
        coutns_banana_means = counts_banana_means / divisor
        counts_juice_stds = counts_juice_stds / divisor
        counts_banana_stds = counts_banana_stds / divisor
    end
    juice_more_juice = counts_juice_means[Int.(counts_juice_means .> counts_banana_means)]
    scatter(counts_juice_means, counts_banana_means, legend=:none, xerr=counts_juice_stds ./ sqrt(106), yerr = counts_banana_stds ./sqrt(106),ylim=(-3.5,3.5),xlim=(-3.5,3.5),lc=:reds,mc=:reds, msc=:reds, guidefontsize=13, tickfontsize=TICKFONTSIZE, legendfontsize=LEGENDFONTSIZE)
    xs = collect(-3:0.1:3)
    plot!(xs, xs,color=:black,linestyle=:dash,alpha=0.7)
    if plot_firing_rate
        xlabel!("Relative firing rate (0.9ml)")
        ylabel!("Relative firing rate (1.5g)")
        title!("Juice vs banana firing rates by neuron")
    else
        xlabel!("Relative spikes juice (0.9ml)")
        ylabel!("Relative spikes banana (1.5g)")
        title!("Juice spikes vs banana spikes by neurons")
    end
    savefig("figures/spikes_scatter_100_400_normalized_3.png")
    # statistics
    float_juice_counts = float.(counts_juice_means)
    float_banana_counts = float.(counts_banana_means)
    corr_coeff = Statistics.cor(float_juice_counts, float_banana_counts)
    p_val = pvalue(CorrelationTest(float_juice_counts, float_banana_counts))
    println("R value $corr_coeff")
    println("p value: $p_val")
    return corr_coeff, p_val
end

function relative_spikes_scatterplot(normalized_counts = true, use_all_situations = false, use_best_fit_line = true,plot_firing_rate=false)
    window_start = 100
    window_end = 400
    window_size = window_end - window_start
    divisor = window_size / 1000 # msec to sec conversion
    println("STARTING FUNCTION")
    if use_all_situations
        println("using all situations")
        counts_juice, counts_banana = neuron_means_std_spikes_in_window_by_situation(neuronlist,[ALL_JUICE_SITUATIONS, ALL_BANANA_SITUATIONS],window_start,window_end,false, normalized_counts, plot_firing_rate)
    else
        println("not doing so")
        counts_juice, counts_banana = neuron_means_std_spikes_in_window_by_situation(neuronlist,[JUICE_SITUATIONS, BANANA_SITUATIONS],window_start,window_end,false, normalized_counts, plot_firing_rate)
    end
    counts_juice_means, counts_juice_stds = counts_juice
    counts_banana_means, counts_banana_stds = counts_banana
    more_juice_juice, more_juice_banana,more_banana_juice,more_banana_banana,more_juice_juice_stds,more_juice_banana_stds,more_banana_juice_stds,more_banana_banana_stds = new_filter_juice_banana(counts_juice_means, counts_banana_means,counts_juice_stds, counts_banana_stds)
    #if plot_firing_rate
    #    more_juice_juice = more_juice_juice / divisor
    #    more_juice_banana = more_juice_banana / divisor
    #    more_banana_juice = more_banana_juice / divisor
    #    more_banana_banana = more_banana_banana / divisor
    #    more_juice_juice_stds = more_juice_juice_stds / divisor
    #    more_juice_banana_stds = more_juice_banana_stds / divisor
    #    more_banana_juice_stds = more_banana_juice_stds / divisor
    #    more_banana_banana_stds = more_banana_banana_stds / divisor
    #end

    scatter(more_juice_juice, more_juice_banana,xerr = more_juice_juice_stds ./ sqrt(106), yerr=more_juice_banana_stds ./sqrt(106),m=:circle,lc=:blue, mc=:blue, msc=:blue,alpha=1,label=false, guidefontsize=11, tickfontsize=TICKFONTSIZE,legendfontsize=LEGENDFONTSIZE, aspect_ratio = :equal)
    scatter!(more_banana_juice, more_banana_banana,xerr = more_banana_juice_stds ./ sqrt(106), yerr=more_banana_banana_stds ./sqrt(106),m=:circle,lc=:orange, mc=:orange, msc=:orange,alpha=1,label=false)
    scatter!(more_juice_juice, more_juice_banana,m=:circle,lc=:black, mc=:black, msc=:black,label=false)
    scatter!(more_banana_juice, more_banana_banana,m=:black,lc=:black, mc=:black, msc=:black,label=false)
    vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
    hline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
    if !plot_firing_rate
        ylims!(-3.5,3.5)
        xlims!(-3.5,3.5)
    else
        ylims!(-5.0, 10.0)
        xlims!(-5.0, 10.0)
    end
    # compute best fit line
    if use_best_fit_line
        if plot_firing_rate
            xs = collect(-5:0.1:10)
        else
            xs = collect(-3:0.1:3)
        end
        float_juice_counts = float.(counts_juice_means)
        float_banana_counts = float.(counts_banana_means)
        corr_coeff = Statistics.cor(float_juice_counts, float_banana_counts)
        p_val = pvalue(CorrelationTest(float_juice_counts, float_banana_counts))
        println("cor coeff: $corr_coeff")
        println("pval: $p_val")
        intercept =1.0
        plot!(xs, (corr_coeff .* xs) .+ intercept,label="Best fit line",color=:red,legend=:topleft)
    end
    if plot_firing_rate
        xlabel!("Relative firing rate juice (0.9ml)")
        ylabel!("Relative firing rate banana (1.5g)")
        title!("Juice vs banana firing rates by neuron")
    else
        xlabel!("Relative spikes juice (0.9ml)")
        ylabel!("Relative spikes banana (1.5g)")
        title!("Juice spikes vs banana spikes by neurons")
    end
    if use_all_situations
        savefig("figures/spikes_scatter_100_400_normalized_colored_correlation_no_best_fit_all_situations_3.png")
    else
        savefig("figures/spikes_scatter_100_400_normalized_colored_correlation_no_best_fit_3.png")
    end
end

function subjective_value_barchart(USE_TITLE=false)
    subjective_values = [0.05, 0.1, 0.5,0.7,1]
    situation_descriptions = ["0.3g Banana", "0.2ml Juice","0.5ml Juice", "1.5g Banana", "0.9ml Juice"]
    bar(situation_descriptions, subjective_values,legend=:none,color=[:orange,:blue,:blue,:orange,:blue],alpha=0.5, guidefontsize=GUIDEFONTSIZE, tickfontsize=12, legendfontsize=LEGENDFONTSIZE)
    xlabel!("Experimental condition")
    ylabel!("Subjective value")
    plot!(size=(750,400))
    if USE_TITLE
        title!("Subjective value by condition")
    end
    savefig("figures/subjective_value_barchart_resized.png")
end

function significance_histogram(use_stored_data = true,use_title = true, place_vline_label = true)
    window_start = 100
    window_end = 400
    N_runs = 1000
    N_significant = 13
    if use_stored_data
        slist = npzread("shuffled_situation_significant_list.npy")
    else
        significant_list = shuffle_by_situation(neuronlist, N_runs,window_start,window_end)
        slist = Int.(significant_list)
        npzwrite("shuffled_situation_significant_list.npy", slist)
    end
    if place_vline_label
        vline_label = "Number of significant coefficients"
    else
        vline_label = ""
    end
    xs = collect(0:N_significant+1)
    xticks = (xs[1:end-1] .+ 0.5, string.(xs))
    histogram(slist,xlim=(0,N_significant + 1),bins=xs,xticks=xticks,bar_width=1,label=false, guidefontsize=GUIDEFONTSIZE+1, tickfontsize=11, legendfontsize=LEGENDFONTSIZE+1,titlefontsize=TITLEFONTSIZE)
    vline!([N_significant + 0.5], linestyle=:dash, alpha=0.9, color=:gray,label=vline_label)
    plot!(size=(500,500))
    xlabel!("Number of significant coefficients")
    ylabel!("Number of samples")
    if use_title
        title!("Empirical Distribution of Significant Coefficients")
    end
    savefig("figures/histogram_shuffled_situation_coefficients_4_resized.png")
end


subjective_value_barchart()
relative_spikes_scatterplot(true, false, true,true)
plot_overlapping_window_timecourses(200,true)
#timecourse_plots_non_overlapping()
all_neurons_condition_histogram(100,400,true,true,true)
significance_histogram(true, false)


# investigating spike distribution and number of trials removed
bucket_count_neuronlist(neuronlist,ALL_SITUATIONS,11000,-3000,8000)
mean_spikes, std_spikes,spikelist = mean_and_std_spikes()

histogram(spikelist)
total_removed = 0
for n in spikelist
    if n < 5
        total_removed +=1
    end
end
total_removed
