# utility functions for plotting for the data analysis
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


normalize(vec) = (vec .- mean(vec) ./ std(vec))

function moving_average(xs, window_size)
    for i in 1:length(xs)
        if i <= window_size
            xs[i] = mean(xs[1:i])
        else
            xs[i] = mean(xs[i-window_size:i])
        end
    end
    return xs
end

function bucket_count_spikes(spikes, window,bucket_start = -2000, bucket_end=3500,use_overlapping_window=false)
    println("In bucket count spikes $use_overlapping_window")
    if use_overlapping_window == true

        return overlapping_window_count_spikes(spikes, window)
    end
    bucket_edges = collect(bucket_start:window:bucket_end)
    counts = zeros(length(bucket_edges))
    for i in 1:length(spikes)
        stimes = spikes[i]
        for time in stimes
            if time <= bucket_edges[1]
                counts[1] +=1
            end
            for b in 1:length(bucket_edges) -2
                if time >= bucket_edges[b] && time <= bucket_edges[b + 1]
                    counts[b+1] +=1
                end
            end
            if time >= bucket_edges[length(bucket_edges)]
                counts[length(bucket_edges)] +=1
            end
        end
    end
    return counts
end

function count_spikes_within_window(spikes, tstart, tend)
    count = 0
    for stimes in spikes    # loop over all spiketrains
        for t in stimes
            if t >= tstart && t<= tend
                count +=1
            end
        end
    end
    return count
end

function overlapping_window_count_spikes(spikes, window_size, bucket_start=-500, bucket_end = 1500,step_size=50,mean_normalize_across_trials = true)
    buckets = []
    bucket_edge = bucket_start
    println("in overlapping window count spikes $bucket_start")
    bucket_starts = collect(bucket_start:step_size:bucket_end)
    counts = zeros(length(bucket_starts))
    for i in 1:length(spikes)
        stimes = spikes[i]
        for time in stimes
            if time <= bucket_starts[1]
                counts[1] +=1
            end
            if time >= bucket_starts[end]
                counts[end] +=1
            end
            for i in 1:length(bucket_starts)
                #if time >= bucket_starts[i] -  && time <= bucket_starts[i] + window_size
                if time >= (bucket_starts[i] - (window_size / 2)) && time <= (bucket_starts[i] + (window_size / 2))
                    counts[i] +=1
                end
            end
        end
    end
    if mean_normalize_across_trials
        L = length(spikes)
        println("mean normalizing across spikes $L")
        counts = counts ./ length(spikes)
    end
    return counts
end


function get_relative_times(spikes, event_times)
    relative_times = []
    for i in 1:length(spikes)
        if typeof(event_times[i]) == Int16
            #println("INDEX: $i")
            push!(relative_times, spikes[i] .- event_times[i])
        else
            e = event_times[i]
            #println("error: $e")
        end
    end
    return relative_times
end

function filter_spikes_events_by_situation(spikes,events,situations, situation_code_list)
    # figure out filter by situatios
    spike_list_situations = []
    event_list_situations = []
    for i in 1:length(spikes)
        if situations[i][2] in situation_code_list
            push!(spike_list_situations, spikes[i])
            push!(event_list_situations, events[i])
        end
    end
    return spike_list_situations, event_list_situations
end

function get_counts_by_situation(spiketimes, events,situation, situation_code, window_size)
    spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes,events, situation, situation_code)
    relative_spikes_filt = get_relative_times(spikes_filt, events_filt)
    counts = bucket_count_spikes(relative_spikes_filt, window_size)
    return counts
end

function get_spikes_per_trial_within_window(spiketimes, events, situations, situation_list,tstart, tend)
    relative_times = get_relative_times(spiketimes, events)
    filtered_situations = []
    trial_counts = zeros(length(relative_times))
    for (i,stimes) in enumerate(relative_times)
        push!(filtered_situations, situations[i])
        for time in stimes
            if time >= tstart && time <= tend
                trial_counts[i] +=1
            end
        end
    end
    return trial_counts, filtered_situations
end

function count_spikes_within_window_by_situation(spiketimes, events, situations, situation_list, tstart, tend)
    spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes,events,situations, situation_list)
    relative_spikes = get_relative_times(spikes_filt, events)
    count = count_spikes_within_window(relative_spikes, tstart, tend)
    return count
end

function bucket_counts_by_situation(spiketimes, events, situations, situation_list, window_size, bucket_start, bucket_end,use_overlapping_window=false)
    spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes, events, situations, situation_list)
    relative_spikes = get_relative_times(spikes_filt, events)
    counts = bucket_count_spikes(relative_spikes, window_size, bucket_start, bucket_end,use_overlapping_window)
    return counts
end


function parse_data_dict(sname="Data/w064-0122.mat")
    data_dict  = matread(sname)["savefile"]["trial"]
    duration = data_dict["duration"]
    spiketimes = data_dict["neuron"]
    L = length(duration)
    bits = data_dict["bit"]
    situations = data_dict["situation"]
    fp_times = []
    stim_onsets = []
    solenoid_1s = []
    solenoid_2s = []
    for i in 1:L
        push!(fp_times, bits[i]["upat"][1])
        push!(stim_onsets, bits[i]["upat"][2])
        push!(solenoid_1s, bits[i]["upat"][5])
        push!(solenoid_2s, bits[i]["upat"][8])
    end
    return spiketimes, bits,situations, fp_times, stim_onsets, solenoid_1s,solenoid_2s,duration
end

function bucket_count_neuronlist(neuronlist, situation_list, window_size, bucket_start, bucket_end)
    all_counts = []
    for neuron in neuronlist
        spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" * ".mat")
        try
            counts = bucket_counts_by_situation(spiketimes, stim_onsets,situations, ALL_SITUATIONS,window_size, bucket_start, bucket_end)
            push!(all_counts, counts)
        catch Exception
            println("ERROR ON NEURON $neuron   $Exception")
        end
    end
    return all_counts
end


function bucket_count_by_situation_neuronlist(neuronlist, situation_list, window_size, bucket_start, bucket_end, crop_end = true,use_overlapping_window=false)
    all_situation_list = []
    for situation in ALL_SITUATIONS
        all_counts = []
        for neuron in neuronlist
            spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" * ".mat")
            try
                counts = bucket_counts_by_situation(spiketimes, stim_onsets,situations, [situation],window_size, bucket_start, bucket_end,use_overlapping_window)
                if crop_end
                    push!(all_counts, counts[2:end-1]) # crop end collapsed buckets
                else
                    push!(all_counts, counts)
                end
            catch e
                println("ERROR ON NEURON $neuron")
                println("Exception $e")
            end
        end
        push!(all_situation_list, all_counts)
    end
    return all_situation_list
end

function count_plot_by_neuron(neuronlist, bucket_start, bucket_end, window_size, normalize_by_prev_window=false,use_overlapping_window=false)
    println("In count plot by neuron $use_overlapping_window")
    all_bucket_counts =bucket_count_by_situation_neuronlist(neuronlist,ALL_SITUATIONS,window_size, bucket_start, bucket_end,true, use_overlapping_window)
    if normalize_by_prev_window == true
        print("not implemented") # I don't think this actually makes sense in this case!
    end
    mean_all_counts = []
    std_all_counts = []
    for i in 1:length(ALL_SITUATIONS)
        counts_matrix = hcat(vcat.(all_bucket_counts[i]...)...)
        mean_counts = mean(counts_matrix, dims=1)
        std_counts = std(counts_matrix, dims=1)
        push!(mean_all_counts, mean_counts)
        push!(std_all_counts, std_counts)
    end
    return mean_all_counts, std_all_counts
end

function filter_situations_by_situationlist(situations, valid_situations)
    situation_list = []
    for situation in situations
        #println(situation)
        if situation[2] in valid_situations
            push!(situation_list, situation)
        end
    end
    return situation_list
end

function filter_situation_by_subjective_value(situations)
    subjective_value_regressor = float.(zeros(length(situations)))
    for i in 1:length(situations)
        sit_val = situations[i][2]
        if sit_val == 5
            subjective_value_regressor[i] = 1
        elseif sit_val == 1
            subjective_value_regressor[i] = 0.7
        elseif sit_val == 4
            subjective_value_regressor[i] = 0.6
        elseif sit_val == 3
            subjective_value_regressor[i] = 0.1
        elseif sit_val == 2
            subjective_value_regressor[i] = 0.05
        end
    end
    return subjective_value_regressor
end

function reward_type_subjective_value_regressor(situations)
    food_regressor = zeros(length(situations))
    juice_regressor = zeros(length(situations))
    for i in 1:length(situations)
        sit_val = situations[i][2]
        if sit_val == 5
            food_regressor[i] = 0.0
            juice_regressor[i] = 1.0
        elseif sit_val == 4
            food_regressor[i] = 0.0
            juice_regressor[i] = 0.6
        elseif sit_val == 3
            food_regressor[i] = 0.0
            juice_regressor[i] = 0.1
        elseif sit_val == 2
            food_regressor[i] = 0.05
            juice_regressor[i] = 0.0
        elseif sit_val == 1
            food_regressor[i] = 0.7
            juice_regressor[i] = 0.0
        end
    end
    return food_regressor, juice_regressor
end


function filter_situations_by_type(situations)
    food_regressor = zeros(length(situations))
    juice_regressor = zeros(length(situations))
    for i in 1:length(situations)
        sit_val = situations[i][2]
        if sit_val in JUICE_SITUATIONS
            juice_regressor[i] = 1
        elseif sit_val in BANANA_SITUATIONS
            food_regressor[i] = 1
        end
    end
    return food_regressor, juice_regressor
end

function situation_bipolar_regressor(situations)
    type_regressor = zeros(length(situations))
    for i in 1:length(situations)
        sit_val = situations[i][2]
        if sit_val in JUICE_SITUATIONS
            type_regressor[i] = 1
        elseif sit_val in BANANA_SITUATIONS
            type_regressor[i] = -1
        end
    end
    return type_regressor
end

function counts_per_trial_window_by_situation(spiketimes, events, situations,situation_list,window_start, window_end, control_normalization = false, per_trial_normalization = false)
    spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes,events,situations, situation_list)
    relative_spikes = get_relative_times(spikes_filt, events_filt)
    # count spikes for each trial
    counts_per_trial = []
    total_control_counts = 0
    Tdiff = window_end - window_start
    for stimes in relative_spikes
        count = 0
        for t in stimes
            if t>= window_start && t<= window_end
                count+=1
            end
            if control_normalization == true
                if t >= window_start - 100 - Tdiff && t <= window_start - 100
                    total_control_counts +=1
                end
            end
        end
        push!(counts_per_trial, count)
    end
    mean_control_counts = total_control_counts / length(relative_spikes)
    if control_normalization == true
        # normalize
        counts_per_trial = counts_per_trial ./ mean_control_counts
    end
    return counts_per_trial
end


function neurons_spike_proportions(neuronlist, situation_1, situation_2, window_start, window_end)
    mean_1s = []
    stds_1s = []
    mean_2s = []
    stds_2s = []
    for neuron in neuronlist
        spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" *".mat")
        spikes_filt_1, events_filt_1 = filter_spikes_events_by_situation(spiketimes,stim_onsets,situations, situation_1)
        relative_spikes_1 = get_relative_times(spikes_filt_1, events_filt_1)
        spikes_filt_2, events_filt_2 = filter_spikes_events_by_situation(spiketimes,stim_onsets,situations, situation_2)
        relative_spikes_2 = get_relative_times(spikes_filt_2, events_filt_2)
        count_1s_per_trial = []
        count_2s_per_trial = []
        Tdiff = window_end - window_start
        for (stimes_1, stimes_2) in zip(relative_spikes_1, relative_spikes_2)
            count_1s = 0
            count_2s = 0
            for (t1,t2) in zip(stimes_1, stimes_2)
                if t1>= window_start && t1<= window_end
                    count_1s+=1
                end
                if t2>= window_start && t2<= window_end
                    count_2s+=1
                end
            end
            if count_1s + count_2s > 0
                # make sure we ignore trials with 0 total spikev (!)
                push!(count_1s_per_trial, count_1s)
                push!(count_2s_per_trial, count_2s)
            end
        end
        totals_per_trial = count_1s_per_trial .+ count_2s_per_trial
        println("TOTALS PER TRIAL: $totals_per_trial")
        proportion_1s_per_trial = count_1s_per_trial ./ totals_per_trial
        println("Proportion 1s: $proportion_1s_per_trial")
        proportion_2s_per_trial = count_2s_per_trial ./ totals_per_trial
        mean_proportion_1s = mean(proportion_1s_per_trial)
        mean_proportion_2s = mean(proportion_2s_per_trial)
        std_proportion_1s = std(proportion_1s_per_trial)
        std_proportion_2s = std(proportion_2s_per_trial)
        push!(mean_1s, mean_proportion_1s)
        push!(mean_2s, mean_proportion_2s)
        push!(stds_1s, std_proportion_1s)
        push!(stds_2s, std_proportion_2s)
    end
    return mean_1s, mean_2s, stds_1s, stds_2s
end




function neuron_means_std_spikes_in_window_by_situation(neuronlist, situation_list,window_start, window_end, control_normalization = false, prev_window_normalize=false, firing_rate_conversion = true)
    counts = []
    for situation in situation_list
        neuron_mean_list =  []
        neuron_std_list = []
        for neuron in neuronlist
            spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" *".mat")
            count_list = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, situation,window_start,window_end,control_normalization)
            count_list = vcat(count_list...)
            mean_count = mean(count_list)
            std_count = std(count_list)
            if prev_window_normalize == true
                prev_window_start = -500
                prev_window_end = 0
                spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" *".mat")
                prev_count_list = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, situation,prev_window_start,prev_window_end,control_normalization)
                prev_count_list = vcat(prev_count_list...)
                prev_mean_count = mean(prev_count_list)
                #println("$mean_count -- $prev_mean_count")
                if firing_rate_conversion
                    divisor = (window_end - window_start) / 1000
                    prev_divisor = (prev_window_end - prev_window_start) / 1000
                    mean_count = mean_count ./ divisor
                    prev_mean_count = prev_mean_count ./ prev_divisor
                end
                mean_count = mean_count.- prev_mean_count
            end
            if firing_rate_conversion
                std_count = std_count ./ divisor
            end
            push!(neuron_mean_list, mean_count)
            push!(neuron_std_list, std_count)
        end
        push!(counts, [neuron_mean_list, neuron_std_list])
    end
    return counts
end


function get_relative_times_2(spikes, event_times,situations)
    relative_times = []
    output_situations = []
    for i in 1:length(spikes)
        if typeof(event_times[i]) == Int16
            push!(relative_times, spikes[i] .- event_times[i])
            push!(output_situations, situations[i])
        else
            e = event_times[i]
        end
    end
    return relative_times, output_situations
end


function subjective_value_reward_type_coeffs(neuronlist,situation_list, start_window = 100, end_window=400,shuffle_situation_flag = false, normalize_regressors = false)
    coeff_list = []
    ts_list = []
    ps_list = []
    #formula = @formula(Spikes ~ Juice1 + Juice2)
    for neuron in neuronlist
        #println("Neuron $neuron")
        spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" *".mat")
        count_list = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, EXPERIMENT_SITUATIONS,start_window,end_window,false)
        filtered_situations = filter_situations_by_situationlist(situations, EXPERIMENT_SITUATIONS)
        spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes,stim_onsets,situations, EXPERIMENT_SITUATIONS)
        relative_spikes, filtered_situations = get_relative_times_2(spikes_filt, events_filt,filtered_situations)
        if shuffle_situation_flag == true
            filtered_situations = shuffle_situations(filtered_situations)
        end
        subjective_value_regressor = filter_situation_by_subjective_value(filtered_situations)
        food_regressor, juice_regressor = filter_situations_by_type(filtered_situations)
        if normalize_regressors == true
            subjective_value_regressor = (subjective_value_regressor .- mean(subjective_value_regressor)) ./ std(subjective_value_regressor)
            food_regressor = (food_regressor .- mean(food_regressor)) ./ std(food_regressor)
            juice_regressor = (juice_regressor .- mean(juice_regressor)) ./ std(juice_regressor)
        end
        data = DataFrame(Spikes=float.(count_list),SubjValue = subjective_value_regressor, Food = food_regressor, Juice=juice_regressor)
        formula = @formula(Spikes ~  SubjValue + Food + Juice)
        linreg_model = lm(formula, data)
        # let's check this fitting actually does soething
        coef_cols = coeftable(linreg_model).cols
        coeffs = coef_cols[1]
        ts = coef_cols[3]
        ps = coef_cols[4]
        push!(coeff_list, coeffs)
        push!(ts_list, ts)
        push!(ps_list, ps)
    end

    return hcat(coeff_list...), hcat(ts_list...), hcat(ps_list...)
end

function reward_type_only_coeffs(neuronlist,situation_list, start_window = 100, end_window=400,shuffle_situation_flag = false, normalize_regressors = false)
    coeff_list = []
    ts_list = []
    ps_list = []
    #formula = @formula(Spikes ~ Juice1 + Juice2)
    for neuron in neuronlist
        #println("Neuron $neuron")
        spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" *".mat")
        count_list = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, EXPERIMENT_SITUATIONS,start_window,end_window,false)
        filtered_situations = filter_situations_by_situationlist(situations, EXPERIMENT_SITUATIONS)
        spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes,stim_onsets,situations, EXPERIMENT_SITUATIONS)
        relative_spikes, filtered_situations = get_relative_times_2(spikes_filt, events_filt,filtered_situations)
        if shuffle_situation_flag == true
            filtered_situations = shuffle_situations(filtered_situations)
        end
        food_regressor, juice_regressor = filter_situations_by_type(filtered_situations)
        if normalize_regressors == true
            food_regressor = (food_regressor .- mean(food_regressor)) ./ std(food_regressor)
            juice_regressor = (juice_regressor .- mean(juice_regressor)) ./ std(juice_regressor)
        end
        data = DataFrame(Spikes=float.(count_list), Food = food_regressor, Juice=juice_regressor)
        formula = @formula(Spikes ~  Food + Juice)
        linreg_model = lm(formula, data)
        # let's check this fitting actually does soething
        coef_cols = coeftable(linreg_model).cols
        coeffs = coef_cols[1]
        ts = coef_cols[3]
        ps = coef_cols[4]
        push!(coeff_list, coeffs)
        push!(ts_list, ts)
        push!(ps_list, ps)
    end

    return hcat(coeff_list...), hcat(ts_list...), hcat(ps_list...)
end

function subjective_value_bipolar_reward(neuronlist,situation_list, start_window, end_window, shuffle_situation_flag = false, normalize_regressors = false)
    coeff_list = []
    ts_list = []
    ps_list = []
    #formula = @formula(Spikes ~ Juice1 + Juice2)
    for neuron in neuronlist
        #println("Neuron $neuron")
        spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" *".mat")
        count_list = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, EXPERIMENT_SITUATIONS,start_window,end_window,false)
        filtered_situations = filter_situations_by_situationlist(situations, EXPERIMENT_SITUATIONS)
        spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes,stim_onsets,situations, EXPERIMENT_SITUATIONS)
        relative_spikes, filtered_situations = get_relative_times_2(spikes_filt, events_filt,filtered_situations)
        if shuffle_situation_flag == true
            filtered_situations = shuffle_situations(filtered_situations)
        end
        subjective_value_regressor = filter_situation_by_subjective_value(filtered_situations)
        type_regressor = situation_bipolar_regressor(filtered_situations)
        if normalize_regressors == true
            subjective_value_regressor = (subjective_value_regressor .- mean(subjective_value_regressor)) ./ std(subjective_value_regressor)
            type_regressor = (type_regressor .- mean(type_regressor)) ./ std(type_regressor)
        end
        data = DataFrame(Spikes=float.(count_list),SubjValue = subjective_value_regressor, Type = type_regressor)
        formula = @formula(Spikes ~  SubjValue + Type)
        linreg_model = lm(formula, data)
        coef_cols = coeftable(linreg_model).cols
        coeffs = coef_cols[1]
        ts = coef_cols[3]
        ps = coef_cols[4]
        push!(coeff_list, coeffs)
        push!(ts_list, ts)
        push!(ps_list, ps)
    end

    return hcat(coeff_list...), hcat(ts_list...), hcat(ps_list...)
end


function reward_type_only_subjective_value_coeffs(neuronlist,situation_list, start_window = 100, end_window=400,shuffle_situation_flag = false, normalize_regressors = false)
    coeff_list = []
    ts_list = []
    ps_list = []
    #formula = @formula(Spikes ~ Juice1 + Juice2)
    for neuron in neuronlist
        #println("Neuron $neuron")
        spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" *".mat")
        count_list = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, EXPERIMENT_SITUATIONS,start_window,end_window,false)
        filtered_situations = filter_situations_by_situationlist(situations, EXPERIMENT_SITUATIONS)
        spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes,stim_onsets,situations, EXPERIMENT_SITUATIONS)
        relative_spikes, filtered_situations = get_relative_times_2(spikes_filt, events_filt,filtered_situations)
        if shuffle_situation_flag == true
            filtered_situations = shuffle_situations(filtered_situations)
        end
        food_regressor, juice_regressor = reward_type_subjective_value_regressor(filtered_situations)
        if normalize_regressors == true
            food_regressor = (food_regressor .- mean(food_regressor)) ./ std(food_regressor)
            juice_regressor = (juice_regressor .- mean(juice_regressor)) ./ std(juice_regressor)
        end
        data = DataFrame(Spikes=float.(count_list), Food = food_regressor, Juice=juice_regressor)
        formula = @formula(Spikes ~  Food + Juice)
        linreg_model = lm(formula, data)
        # let's check this fitting actually does soething
        coef_cols = coeftable(linreg_model).cols
        coeffs = coef_cols[1]
        ts = coef_cols[3]
        ps = coef_cols[4]
        push!(coeff_list, coeffs)
        push!(ts_list, ts)
        push!(ps_list, ps)
    end

    return hcat(coeff_list...), hcat(ts_list...), hcat(ps_list...)
end


function shuffle_situations(situations)
    shuffled_list = []
    for situation in situations
        situation[2] = Int(trunc(rand(Uniform(1,6))))
        push!(shuffled_list, situation)
    end
    return shuffled_list
end

function N_significant(ps,p_idx, p_value)
    ps_vec = ps[p_idx, :]
    N = 0
    for p in ps_vec
        if p <= p_value
            N+=1
        end
    end
    return N
end


function shuffle_by_situation(neuronlist, N_shuffles, window_start, window_end)
    N_significants = []
    for i in 1:N_shuffles
        coeff_list, ts_list, ps_list = subjective_value_bipolar_reward(neuronlist, EXPERIMENT_SITUATIONS, window_start, window_end, true)
        N = N_significant(ps_list, 3, 0.05)
        push!(N_significants, N)
    end
    return N_significants
end

function window_size_shuffle(neuronlist, N_shuffles, avoid_start, avoid_end, min_val, max_val)
    tdiff = avoid_end - avoid_start
    N_significants = []
    for i in 1:N_shuffles
        start = Int(trunc(rand(Uniform(min_val, max_val))))
        while ((start > avoid_start) && (start < avoid_end))
            start = Int(trunc(rand(Uniform(min_val, max_val))))
        end
        println("start : $start")
        window_end = start + tdiff
        coeff_list, ts_list, ps_list = subjective_value_bipolar_reward(neuronlist,EXPERIMENT_SITUATIONS,start, window_end)
        N = N_significant(ps_list, 3, 0.05)
        push!(N_significants, N)
    end
    return N_significants
end

function count_greater_than(l, geval)
    N = 0
    for el in l
        if el >= geval
            N+=1
        end
    end
    return N
end


function new_filter_juice_banana(counts_juice_means, counts_banana_means,counts_juice_stds, counts_banana_stds)
    println("counts juice means: $counts_juice_means")
    more_juice_juice = []
    more_juice_banana = []
    more_banana_juice = []
    more_banana_banana = []
    more_juice_juice_stds = []
    more_juice_banana_stds = []
    more_banana_juice_stds = []
    more_banana_banana_stds = []
    for i in 1:length(counts_juice_means)
        count_juice = counts_juice_means[i]
        count_banana = counts_banana_means[i]
        if count_juice > count_banana
            push!(more_juice_juice, count_juice)
            push!(more_juice_banana, count_banana)
            push!(more_juice_juice_stds, counts_juice_stds[i])
            push!(more_juice_banana_stds, counts_banana_stds[i])
        else
            push!(more_banana_juice, count_juice)
            push!(more_banana_banana, count_banana)
            push!(more_banana_juice_stds, counts_juice_stds[i])
            push!(more_banana_banana_stds, counts_banana_stds[i])
        end
    end
    return more_juice_juice, more_juice_banana,more_banana_juice,more_banana_banana,more_juice_juice_stds,more_juice_banana_stds,more_banana_juice_stds,more_banana_banana_stds
end
