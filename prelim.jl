## reanalysis  to figure out what's going on
# once again julia is best for this analysis
# and I can base analysis on prev code
# all I basically have to do is go in figure 4, figure out which situations are which
# it's just plotting some kind of binned version of them
# from each money -- and assuming that they are the same across trials
# not a given and then from this you just have to see if there are differences between neurons
# while they averaged across neurons right

# right ARGH I don't know what is wrong with me. i'm just slowly dying and getting crushed by STUPID
# executive functionoing issues and thinking too much about the future
# but being unable to actually do / procrastinating awfully out of anything to do in the present
# which is working so damn badly. I'm MEANT TO be able to do ruthless short term execution in favour of a long term objective
#but at the moment I completely lack both a short and a long term objective in a way which makes these things easy
# which is super sucking so badly. I can't even get executive functioning things that I NEED to do done
# which is really really hurting me massively and I want to get all of my papers out ASAP
# and I haven't got anything out or written anything or done anything for MONTHS
# I was struggling with this so bad aroudn neurips time and did okay there briefly but only because I had already done anything
# my work this month is jsut SO SO OS BAD this year and has been c ompletely dominated with a.) getting all sorts of stuff done and out
# and b.) just failing miserably to get any papers done and then ot move on etc
# it's really hurting me so so so badly. I'm too spread out. I'm doing to o many things. I need to a.) get a move on with executive functioning
# since I'm just SO SO SO BAD AT THIS. I NEED TO GET PAID BY VERSES and I NEED TO DO THE VERSES THINGS
# this is just ridicul0ous I need to sort out the tax things, sort out all the vairous things things
# and sort out investments etc is also obvous
# but above all I just need to get out these paper changes. These revisions are KILLER to me because they are so fucking boring and I don't want to do them
# and I'm not excited at all and I just procrastinate but then it kills me from doing any proper work on them.
# like tomorrow, pretty much  right GOT to get through everything and it's so fucking difficult...
# hopefully tomorrow don't have any meetings either so just need to STOP PROCRASTINATING, SIT DOWN AND FIX ALL OF THESE STUPID ISSUES
# ONE AFTER THE OTHER. And will hopefully DO THIS ONE TONIGHT
# which is necessary at least but FUCK. Already wasted ANOTHER HOUR and NOT in a good situation with this so going to be very very stuck
# and super duper stressed out SHOULD HAVE DONE THIS YESTERDAY -- I GUESS I HAVE 1.5 hours or so to actually make this work
# as going to be working VERY late tonight for deadline and that's not even able to do Yuhangs shit. Like this is the issue
# I'm spending SO LONG on speculative stuff and just flailing becaues I'm getting crushed under such stress
# that it's just really killing me that I'm not doing anyhting even though I'm getting curhsed and I have IWAI deadlines in 10 days (!!!) so need to haev ab unch of papers written
# and need to navigate all this craziness in my personal life and it's just killing me ARGH
# and daily 10am meetings with Yuhang are killing me and not using the weekend and having TOO MANY PAPERS and EVEYRTHIGN IS KILLING ME I DON@T KNOW WHAT TO DO AND IM NOT INTERESTED IN IT ARGH
# right and I hate being in here and I'm just SO SO SO STRESSSED
# let's just get this fucking thing done. It's not even that much it's just a bunch of itty bitty issues
# I can do it all in a DAY if I just FOCUS and actually PUT MY HEAD DOWN and just FUCKING DO IT. This is actually ridiculous and I don't know why I procrastinate so much argh


# right let's refresh my memory of where i got to with this and what I need to do
# basically I need to a.) run a linear regression to understand the coefficients
# and try to remember what the linreg was. b.) I need to do the weird shuffling thing and prove its significant
# if I can do this then that would be great for the paper and then I can just try to finish off RB paper as well
# GOD I hate Rafal and his endless ifddling with papers and I WISH that I could just get all the executive functioning stuff DONE
# I just need a basic TODO list every day. That's literally IT and just do stuff off there I don't even have basic stuff and I'm killing myself
# by not having any systems and just not being into this stuff anymore. I FUCKING HATE THIS PAPER
# and I hate lal the rafal papers and they are done thye have been done for over a year and it's just too damn long
# I'm not having it it's just insane and so fucking boring and I hate it. let's just get it done!
using MAT
using Plots
using GLM
using DataFrames
using Distributions
using StatsPlots

# right so let's at least get the additional things here
# see if I can take it easy with the work and don't do any work you don't have to do
# I'll take this meeting from in bed I think so we'll see waht's going on


## TODO: ask Rafal/Armin if I can share this data when paper is released


cd("Desktop/work/rb_data_analysis/second_armin_analysis")

data_dict = matread("Data/w065-0373.mat")["savefile"]["trial"]

duration = data_dict["duration"]
spiketimes = data_dict["neuron"]
L = length(duration)
bits = data_dict["bit"]
fp_times = []
stim_onsets = []
for i in 1:L
    push!(fp_times, bits[i]["upat"][1])
    push!(stim_onsets, bits[i]["upat"][2])
end

fp_times

for i in 1:10
    xs = [i for i in 1:length(spiketimes[i])]
    plot(xs,spiketimes[i][:])
end
xs = [i for i in 1:length(spiketimes[1])]
spiketimes[1][:]
plot(xs,spiketimes[1][:])
maximum(duration)


# let's have a look at how this functions
relative_times = []
for i in 1:L
    if typeof(fp_times[i]) == Int16
        push!(relative_times, spiketimes[i] .- fp_times[i])
    end
end
relative_times
# this is looking much more reasonable -- so let's have al ook

bucket_edges = [-1000,-500,0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000]
counts = zeros(length(bucket_edges))
for i in 1:length(relative_times)
    stimes = relative_times[i]
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
            counts[end] +=1
        end
    end
end

xs = [b for b in bucket_edges]
bar(xs, counts)
# hyh fair enough -- so quite a lot of spikes throughout actually
# which is odd. Idk if the trial or something has changed
# need to actually look at and analyse this
# so why is there a large count right at the end

stim_relative_times = []
for i in 1:L
    if typeof(stim_onsets[i]) == Int16
        push!(stim_relative_times, spiketimes[i] .- stim_onsets[i])
    end
end
stim_relative_times


# it's very nice already having code to do all of this
# I just need to figure out what is what
situations = data_dict["situation"]


function overlapping_window_count_spikes(spikes, window_size, bucket_start=-500, bucket_end = 1500,step_size=10)
    buckets = []
    bucket_edge = bucket_start
    println("in overlapping window count spikes $bucket_start , $bucket_end, $step_size")
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
                if time >= bucket_starts[i] && time <= bucket_starts[i] + window_size
                    counts[i] +=1
                end
            end
        end
    end
    return counts
end


function bucket_count_spikes(spikes, window,bucket_start, bucket_end,use_overlapping_window=false)
    println("In bucket count spikes $use_overlapping_window")
    if use_overlapping_window == true

        return overlapping_window_count_spikes(spikes, window,bucket_start, bucket_end)
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
bucket_start = -2500
bucket_end = 8000
spikes = bucket_count_spikes(stim_relative_times, 100,bucket_start, bucket_end)
xs = collect(bucket_start:100:bucket_end)
plot(xs, spikes)

# relative times function
function get_relative_times(spikes, event_times)
    relative_times = []
    #println("HELLO")
    for i in 1:length(spikes)
        if typeof(event_times[i]) == Int16
            #println("INDEX: $i")
            push!(relative_times, spikes[i] .- event_times[i])
        else
            e = event_times[i]
            #println("argh: $e")
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

normalize(vec) = (vec .- mean(vec) ./ std(vec))

ALL_SITUATIONS = [1,2,3,4,5,25]
spike_counts, filtered_situations = get_spikes_per_trial_within_window(spiketimes, stim_onsets,situations, [2],100,1000)
spike_counts = normalize(spike_counts)
plot(spike_counts)
# huh so these are across trials!? I just need to figure out what to do here
# but making some progress
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

# this is the function I want -- then just mean and std
# to get this solved!
counts = bucket_counts_by_situation(spiketimes, stim_onsets, situations, [5], 100, bucket_start, bucket_end)
plot(counts)
#count_spikes_within_window_by_situation(spiketimes,stim_onsets, situations, ALL_SITUATIONS,100,1000)


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

# I think it would be good at this point since it's looking close to completion
# to actually organize this code and make it ince
# that would be the goal

spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w064-0124.mat")

neuronlist = [359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377]

function bucket_count_neuronlist(neuronlist, situation_list, window_size, bucket_start, bucket_end)
    all_counts = []
    for neuron in neuronlist
        spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0$neuron" * ".mat")
        try
            counts = bucket_counts_by_situation(spiketimes, stim_onsets,situations, ALL_SITUATIONS,window_size, bucket_start, bucket_end)
            push!(all_counts, counts)
        catch Exception
            println("ERROR ON NEURON $neuron")
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
            catch Exception
                println("ERROR ON NEURON $neuron")
            end
        end
        push!(all_situation_list, all_counts)
    end
    return all_situation_list
end


juice_neuron = [359]
food_neuron = [368]
bucket_start = -500
bucket_end = 1500
window_size = 100
standard_neuron = [366]
# 364 and 366 seem reasonable
# I don't think normalizing by the previous window is very useful HERE but probably will be for the earlier scatter plot
# if I can get this sorted then watch another physics lect it will mean no physics tomorrow (or sun I guess) and so can make a lot of progress when not drunk/celebrating
# WHY DO THESE GET SO LARGE!!? SHOULDNT HAVE 25 SPIKES
function count_plot_by_neuron(neuronlist, bucket_start, bucket_end, window_size, normalize_by_prev_window=false,use_overlapping_window=false)
    println("In count plot by neuron $use_overlapping_window")
    all_bucket_counts =bucket_count_by_situation_neuronlist(neuronlist,ALL_SITUATIONS,window_size, bucket_start, bucket_end,true, use_overlapping_window)
    if normalize_by_prev_window == true
        print("not implemented") # yeah I don't think this actually makes sense in this case!
    end
    mean_all_counts = []
    std_all_counts = []
    for i in 1:length(ALL_SITUATIONS)
        counts_matrix = hcat(vcat.(all_bucket_counts[i]...)...)
        # TODO
        mean_counts = mean(counts_matrix, dims=1)
        std_counts = std(counts_matrix, dims=1)
        push!(mean_all_counts, mean_counts)
        push!(std_all_counts, std_counts)
    end
    return mean_all_counts, std_all_counts
end
# I don't know why I'm struggling to focus so much

using Statistics

function moving_average(xs, window_size)
    for i in 1:length(xs)
        if i <= window_size
            xs[i] = mean(xs[1:i])
        #elseif i > length(xs) - window_size
        #    xs[i] = mean(xs[i:end])
        else
            print("in middle thing: $i")
            xs[i] = mean(xs[i-window_size:i])
        end
    end
    return xs
end


plot(moving_average(xs,4))
bucket_start = -500
bucket_end = 1000
window_size = 200
standard_neuron = [366]
xs = collect(bucket_start:10:bucket_end)
juice_means, juice_stds = count_plot_by_neuron(juice_neuron,bucket_start, bucket_end, window_size,true,true)
food_means, food_stds = count_plot_by_neuron(food_neuron, bucket_start, bucket_end, window_size,true,true)
standard_means, standard_stds = count_plot_by_neuron(standard_neuron, bucket_start, bucket_end, window_size,true,true)
plot(xs[2:end-1],juice_means[5][1,:], label="0.9ml juice",color=:blue,tickfontsize=10)
plot!(xs[2:end-1],juice_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
plot!(xs[2:end-1],juice_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
plot!(xs[2:end-1],juice_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
plot!(xs[2:end-1],juice_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
xlabel!("Time (ms) after cue")
ylabel!("Spikes within window")

title!("Spikes of Juice Responsive Neuron")
savefig("juice_neuron_all_conditions_overlapping_window_3.png")

juice_means

plot(xs[2:end-1],food_means[5][1,:], label="0.9ml juice",color=:blue)
plot!(xs[2:end-1],food_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
plot!(xs[2:end-1],food_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
plot!(xs[2:end-1],food_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
plot!(xs[2:end-1],food_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
xlabel!("Time (ms) after cue")
ylabel!("Spikes within window")
title!("Spikes of Food Responsive Neuron")
savefig("food_neuron_all_conditions_overlapping_window_plot_2.png")
# can't do a ribbon since all stds are
plot(xs[2:end-1],standard_means[5][1,:], label="0.9ml juice",color=:blue)
plot!(xs[2:end-1],standard_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
plot!(xs[2:end-1],standard_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
plot!(xs[2:end-1],standard_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
plot!(xs[2:end-1],standard_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
xlabel!("Time (ms) after cue")
ylabel!("Spikes within window")
title!("Spikes of Value Responsive Neuron")
savefig("value_only_neuron_all_conditions_window_plot.jpg.png")


xs = collect(bucket_start:window_size:bucket_end)
juice_means, juice_stds = count_plot_by_neuron(juice_neuron,bucket_start, bucket_end, window_size,true)
food_means, food_stds = count_plot_by_neuron(food_neuron, bucket_start, bucket_end, window_size,true)
standard_means, standard_stds = count_plot_by_neuron(standard_neuron, bucket_start, bucket_end, window_size,true)
juice_stds
mean_all_counts
# yeah I'm not sure I can make a nice looking graph of this
# how should I do this?
#5 = 0.9ml juice
#1 = 1.5kg bananas
# 4 = 0.5ml juice
# 3 = 0.2ml juice
# 2  = 0.3kg bananas
# okay well they don't look super nice but at least some of them i.e. the most dramatic ones can be one
# I'm actually getting increasinfly convinced that there is an effect here.  -- will need to do some messing aroudn with the model
# to see what can be done to improve this but at least done the RB stuff and mostly the UHN stuff but got LOADS more to do tomorrow
# basically the temporal PC only

# yeah so this is the issue right, these are the NEURONS in the SAME condition
# so it's misleading to do it across plots which RAFAL thinks it is
# perhaps a better plot is response of juice neuron across plots
# okay this moving average is not working -- WHY cannot you FOUCS

plot(juice_means[5][1,:])
plot!(moving_average(juice_means[5][1,:],5))
# okay this completely distorts the pattern
xs[2:end-1]
juice_means[5][1,:]
plot(xs[2:end-1],juice_means[5][1,:], label="0.9ml juice",color=:blue)
plot!(xs[2:end-1],juice_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
plot!(xs[2:end-1],juice_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
plot!(xs[2:end-1],juice_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
plot!(xs[2:end-1],juice_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
xlabel!("Time (ms) after cue")
ylabel!("Spikes within window")
title!("Spikes of Juice Responsive Neuron")
savefig("juice_neuron_all_conditions_window_plot_2.jpg")

plot(xs[2:end-1],food_means[5][1,:], label="0.9ml juice",color=:blue)
plot!(xs[2:end-1],food_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
plot!(xs[2:end-1],food_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
plot!(xs[2:end-1],food_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
plot!(xs[2:end-1],food_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
xlabel!("Time (ms) after cue")
ylabel!("Spikes within window")
title!("Spikes of Food Responsive Neuron")
savefig("food_neuron_all_conditions_window_plot_2.jpg")

plot(xs[2:end-1],standard_means[5][1,:], label="0.9ml juice",color=:blue)
plot!(xs[2:end-1],standard_means[4][1,:], label="0.5ml juice",color=:blue,linestyle=:dashdot)
plot!(xs[2:end-1],standard_means[3][1,:], label="0.3ml juice",color=:blue,linestyle=:dot)
plot!(xs[2:end-1],standard_means[1][1,:], label="1.5g banana",color=:orange,linestyle=:solid)
plot!(xs[2:end-1],standard_means[2][1,:], label="0.3g banana",color=:orange,linestyle=:dashdot)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label="Stimulus Onset")
xlabel!("Time (ms) after cue")
ylabel!("Spikes within window")
title!("Spikes of Value Responsive Neuron")
savefig("value_only_neuron_all_conditions_window_plot.jpg")
# so there's still these same kinds of issues and very strong value coding essentially
# so this is vaguely encouraging

plot(xs[2:end-1],juice_means[1][1,:], label="Juice Neuron",color=:blue)
plot!(xs[2:end-1], food_means[1][1,:], label="Food Neuron",color=:yellow)
plot!(xs[2:end-1], juice_means[1][])
xlabel!("Time (ms) after cue")
ylabel!("Spikes within window")
title!("Spikes after receiving 1.5g Bananas")
savefig("food_neuron_window_plot.jpg")

##okay I'm working on this let's get Rafal's additional plots done
# now that I've FINALLY done UHN camera ready


plot(xs[2:end-1],juice_means[5][1,:], label="Juice Neuron")
plot!(xs[2:end-1], food_means[5][1,:], label="Food Neuron")
xlabel!("Time (ms) after cue")
ylabel!("Spikes within window")
title!("Spikes after receiving 0.9ml Juice")
savefig("juice_neuron_window_plot.jpg")

plot!(xs[2:end-1],juice_means[2][1,:], label=ALL_SITUATIONS[2])
plot!(xs[2:end-1],juice_means[3][1,:], label=ALL_SITUATIONS[3])
plot!(xs[2:end-1],juice_means[4][1,:], label=ALL_SITUATIONS[4])
plot!(xs[2:end-1],juice_means[5][1,:], label=ALL_SITUATIONS[5])
#plot!(xs,mean_all_counts[6][1,:], label=ALL_SITUATIONS[6])

# so this is nice, we actually see a juice vs food neuron situation occur
# it's quite a bit noisier than his one  but there are clear peaks
# let's think about how to address this one properly
# because that way I can have all results right there for me
# right let's reorder the situations
# right this is a good stopping point -- can recreate all the figures
# now let's just continue tomorrow!

# situations
#5 = 0.9ml juice subj value 1
#1 = 1.5kg bananas 0.7
# 4 = 0.5ml juice 0.6
# 3 = 0.2ml juice 0.1
# 2  = 0.3kg bananas 0.05

# let's construct the bar chart for this
USE_SUBJ_VALUE_TITLE = false
subjective_values = [0.05, 0.1, 0.5,0.7,1]
situation_descriptions = ["0.3g Banana", "0.2ml Juice","0.5ml Juice", "1.5g Banana", "0.9ml Juice"]
bar(situation_descriptions, subjective_values,legend=:none,tickfontsize=10,color=[:orange,:blue,:blue,:orange,:blue],alpha=0.5)
xlabel!("Experimental Condition")
ylabel!("Subjective Value")
if USE_SUBJ_VALUE_TITLE
    title!("Estimated Subjective Value")
end
savefig("subjective_value_barchart.png")

# juice situations = 5,4,3
# 1,2 = banana situations
# so what

# right so I need a scatter plot and then just generally to fix this up
ALL_JUICE_SITUATIONS = [3,4,5]
ALL_BANANA_SITUATIONS = [1,2]
JUICE_SITUATIONS = [5]
BANANA_SITUATIONS = [1]


spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0359" *".mat")
# right let's create our filter spikes by situations
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
    # this one is a bit bullshit since it's read off by eye from an incorrect chart
    # but that's just the way of it and it'll probably screw up my regression
    # I'm not sure what I should do to do this properly
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

ALL_JUICE_SITUATIONS = [3,4,5]
ALL_BANANA_SITUATIONS = [1,2]



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


# let's think about htis. What do I want? I basically want, for each situation, the proportion of total spikes in the juice vs the food condition
# then if there are like 70% spikes in juice or something things are looking good
# let's try to get this down

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




function neuron_means_std_spikes_in_window_by_situation(neuronlist, situation_list,window_start, window_end, control_normalization = false, prev_window_normalize=false)
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
                println("$mean_count -- $prev_mean_count")
                mean_count = mean_count.- prev_mean_count
            end
            push!(neuron_mean_list, mean_count)
            push!(neuron_std_list, std_count)
        end
        push!(counts, [neuron_mean_list, neuron_std_list])
    end
    return counts
end

bar_xticks = (collect(1:L),string.(neuronlist)[1,:] )
bar(mean_1s,yerr = stds_1s ./ L, legend = :none, xticks = bar_xticks)
xlabel!("Neuron")
ylabel!("Spike Proportion Juice")
title!("Proportion of total spikes to juice or food reward")
savefig("figures/proportions_bar_100_400.jpg")

bar_xticks = (collect(1:L),string.(neuronlist)[1,:] )
bar(mean_1s .- mean(mean_1s),yerr = stds_1s ./ L, legend = :none, xticks = bar_xticks)
xlabel!("Neuron")
ylabel!("Proportion Difference from mean response")
title!("Proportion Differences to juice or food reward")
savefig("figures/proportions_diff_100_400.jpg")

xs = [1,2,3]
xticks = ["1", "2" ,"3"]
bar(xs, xticks=xticks)

# right let's figure this out
EXPERIMENT_SITUATIONS = [1,2,3,4,5]


function get_relative_times_2(spikes, event_times,situations)
    relative_times = []
    output_situations = []
    #println("HELLO")
    for i in 1:length(spikes)
        if typeof(event_times[i]) == Int16
            #println("INDEX: $i")
            push!(relative_times, spikes[i] .- event_times[i])
            push!(output_situations, situations[i])
        else
            e = event_times[i]
            #println("argh: $e")
        end
    end
    return relative_times, output_situations
end



spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0360" *".mat")
count_list = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, EXPERIMENT_SITUATIONS,100,400,false)
filtered_situations = filter_situations_by_situationlist(situations, EXPERIMENT_SITUATIONS)
spikes_filt, events_filt = filter_spikes_events_by_situation(spiketimes,stim_onsets,situations, EXPERIMENT_SITUATIONS)
events_filt
relative_spikes, output_situations = get_relative_times_2(spikes_filt, events_filt, filtered_situations)
subjective_value_regressor = filter_situation_by_subjective_value(filtered_situations)
food_regressor, juice_regressor = filter_situations_by_type(filtered_situations)
linreg_model = lm(formula, data)
# idk why this has never been problematic before

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

coeff_list, ts_list, ps_list= subjective_value_reward_type_coeffs(neuronlist, ALL_SITUATIONS, 100,400, false, true)
npzwrite("banana_juice_reward_type_coeff_list_normalized.npy", coeff_list)
npzwrite("banana_juice_reward_type_ts_list_normalized.npy", ts_list)
npzwrite("banana_juice_reward_type_ps_list_normalized.npy", ps_list)


### plot scatter plot of the coeffs
juice_coeffs = coeff_list[3,:]
food_coeffs = coeff_list[4,:]
scatter(juice_coeffs, food_coeffs, legend=:none ,ylim=(-1,2),xlim=(-1,2))

xs = collect(-2:0.1:2)
#plot!(xs, xs,color=:black,linestyle=:dash,alpha=0.7)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
hline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
xlabel!("Juice Coefficient")
ylabel!("Food Coefficient")
title!("Food vs Juice Coefficients")
savefig("figures/food_juice_coefficients_scatter.png")


coeff_list, ts_list, ps_list = subjective_value_bipolar_reward(neuronlist, ALL_SITUATIONS,100,400,false, true)
npzwrite("banana_juice_bipolar_coeff_list_normalized.npy", coeff_list)
npzwrite("banana_juice_bipolar_ts_list_normalized.npy", ts_list)
npzwrite("banana_juice_bipolar_ps_list_normalized.npy", ps_list)

using NPZ
coeff_list, ts_list, ps_list = reward_type_only_coeffs(neuronlist, ALL_SITUATIONS,100,400,false, true)
npzwrite("reward_type_only_coeff_list_normalized.npy", coeff_list)
npzwrite("reward_type_only_ts_list_normalized.npy", ts_list)
npzwrite("reward_type_only_ps_list_normalized.npy", ps_list)

juice_coeffs = coeff_list[2,:]
food_coeffs = coeff_list[3,:]
scatter(juice_coeffs, food_coeffs, legend=:none ,ylim=(-1,2),xlim=(-1,2))

xs = collect(-2:0.1:2)
#plot!(xs, xs,color=:black,linestyle=:dash,alpha=0.7)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
hline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
xlabel!("Juice Coefficient")
ylabel!("Food Coefficient")
title!("Food vs Juice Coefficients")
savefig("figures/reward_type_only_food_juice_coefficients_scatter.png")



### reward type subjective vcalue coeffs
coeff_list, ts_list, ps_list = reward_type_only_subjective_value_coeffs(neuronlist, ALL_SITUATIONS,100,400,false, true)
npzwrite("reward_type_only_subjective_value_coeff_list_normalized.npy", coeff_list)
npzwrite("reward_type_only_subjective_value_ts_list_normalized.npy", ts_list)
npzwrite("reward_type_only_subjective_value_ps_list_normalized.npy", ps_list)

juice_coeffs = coeff_list[2,:]
food_coeffs = coeff_list[3,:]
scatter(juice_coeffs, food_coeffs, legend=:none ,ylim=(-1,2),xlim=(-1,2))

xs = collect(-2:0.1:2)
#plot!(xs, xs,color=:black,linestyle=:dash,alpha=0.7)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
hline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
xlabel!("Juice Coefficient")
ylabel!("Food Coefficient")
title!("Food vs Juice Coefficients")
savefig("figures/reward_type_only_subjective_value_food_juice_coefficients_scatter.png")

# okay, so this is encouraigng, it looks like we slightly improve the size of the value regressors without messing with the significance too much which sounds about right
# although they also have strong effects of value which is very interesting regardless.
# okay this is good. so I have these plots. what was the other thing he wanted messing with
# right, so I was thinking about this and need to plot how many coeffs are significant
function N_significant(ps,p_idx, p_value)
    ps_vec = ps[p_idx, :]
    #println(ps_vec)
    N = 0
    for p in ps_vec
        #println(p)
        if p <= p_value
            N+=1
        end
    end
    return N
end

ps_list[2,:]
N_significant(ps_list, 3, 0.05)

# right let's see if it works. I bet it'll be completely average and we'll have to avoid this
# but we will see... not totally sure as to proper thing to use

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

# let's see if this can do anything
# will get some nice results in here at the
using NPZ
USE_EMPIRICAL_COEFFS_TITLE = false
significant_list = shuffle_by_situation(neuronlist, 1000,100,400)
slist = Int.(significant_list)
npzwrite("shuffled_situation_significant_list.npy", slist)
slist = npzread("shuffled_situation_significant_list.npy")
#xs = collect(0:11)
xs = collect(0:13)
xticks = (xs[1:end-1] .+ 0.5, string.(xs))
histogram(slist,xlim=(0,12),bins=xs,xticks=xticks,bar_width=1,label=false,guidefontsize=15,legendfontsize=10,tickfontsize=10,titlefontsize=15)
vline!([10.5], linestyle=:dash, alpha=0.9, color=:gray,label="p < 0.01")
xlabel!("Number of Significant Coefficients")
ylabel!("Number of Samples")
if USE_EMPIRICAL_COEFFS_TITLE
    title!("Empirical Distribution of Significant Coefficients")
end
savefig("histogram_shuffled_situation_coefficients_4.jpg")
# oh wow this is really working well it seems with random situations
# fantastic, this actually appears to be working -- will build up the histogram of this right now
# seems reasonable
# okay so this is working well. I can nearly construct a figure with this data analysis stuff
# the PC theory stuff isn't so bad just a bit of figure fiddling -- I might have to recompute some things
#which is highly irritating but not the end of the world... essentially you should ALWAYS SAVE YOUR DATA SO YOU DONT HAVE TO RUN AMSSIVE COMPUTATIONS
# OVER AND OVER AGAIN TO CHANGE PLOT WORDING -- SO STUPID!

significant_list = window_size_shuffle(neuronlist, 1000, 100,400,-1000,8000)
xs = collect(1:11)
npzwrite("significant_list_3.npy", Int.(significant_list))
slist = Int.(significant_list)
histogram(xs,slist,legend=:none,bins=xs)
xlabel!("Number of Significant Coefficients")
ylabel!("Number of Samples")
title!("Empirical Distribution of Significant Coefficients")
savefig("histogram_empirical_coefficients.jpg")
count_greater_than(slist, 10)
maximum(slist)

# okay so this histogram is doing well shuffling over situation, no major issues
# right, I think whattitle!("Empirical Distribution of Significant Coefficients")
 #rafal is wanting is a histogram over counts for a specific (highly angry) neuron

# let's take neuron 359 as our heavy juice neuron and neuron 368 as our heavy juice neuron
# and let's plot the distribution of trials of each in the high juice or high food condition
# histogramming them. This SHOULD be doable. I really like being able to do interactive statistics and visualize liek this
# it's like a better julia notebook and makes data analysis easy if absurdly messy!

# it's actually hilarious how so much time goes into making plotting appear decent and figuring out how to plot overlapping histograms in julia
# I mean I'm also highly confused by pyton. I can just never get roudn to understandinf it!

spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0359" *".mat")
count_list_juice_neuron_juice = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, JUICE_SITUATIONS,100,400,false)
count_list_juice_neuron_food = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, BANANA_SITUATIONS, 100,400,false)
xs = collect(1:11)
count_list_juice_neuron_juice
histogram(xs,count_list_juice_neuron_juice,nbins=10,alpha=0.5,label="Juice Response")
histogram!(xs,count_list_juice_neuron_food,nbins=10, alpha=0.5, label="Food Response")
xlabel!("Spike Count in Window")
ylabel!("Number of Trials")
title!("Juice Neuron Response to Food/Juice")
savefig("juice_neuron_response_histogram.jpg")

#5 = 0.9ml juice
#1 = 1.5kg bananas
# 4 = 0.5ml juice
# 3 = 0.2ml juice
# 2  = 0.3kg bananas

spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0359" *".mat")
situation_1 = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, [1],100,400,false)
situation_2 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [2], 100,400,false)
situation_3 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [3], 100,400,false)
situation_4 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [4], 100,400,false)
situation_5 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [5], 100,400,false)
histogram(xs,[situation_1,situation_2,situation_3,situation_4,situation_5],nbins=10,alpha=0.5,label=["1.5g Banana" "0.3g Banana" "0.2ml Juice" "0.5ml Juice" "0.9ml Juice"],layout=(5,1),xlims=(0,9),ylims=(0,6),xticks=xs)


using Plots.PlotMeasures
# FOOD NEURON
spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0368" *".mat")
situation_2_1 = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, [1],100,400,false)
situation_2_2 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [2], 100,400,false)
situation_2_3 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [3], 100,400,false)
situation_2_4 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [4], 100,400,false)
situation_2_5 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [5], 100,400,false)

# STANDARD NEURON
spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0366" *".mat")
situation_3_1 = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, [1],100,400,false)
situation_3_2 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [2], 100,400,false)
situation_3_3 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [3], 100,400,false)
situation_3_4 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [4], 100,400,false)
situation_3_5 = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, [5], 100,400,false)

xs = collect(0:10)

xs = collect(0:10)
xticks = (xs[1:end-1] .+ 0.5, string.(xs))

#hist_xlabel = "Spike Count in Window"
#hist_ylabel = "Number of Trials"
hist_xlabel = ""
hist_ylabel = ""
hist_plot = histogram(xs,[situation_1,situation_2_1,situation_2, situation_2_2,situation_3,situation_2_3,situation_4,situation_2_4,situation_5,situation_2_5],alpha=0.5,label=["1.5g Banana" "1.5g Banana" "0.3g Banana" "0.3g Banana" "0.2ml Juice" "0.2ml Juice" "0.5ml Juice" "0.5ml Juice" "0.9ml Juice" "0.9ml Juice"],layout=(5,2),xlims=(0,10),ylims=(0,6),xticks=xticks,size=(1000,800),xlabel=hist_xlabel, ylabel=hist_ylabel, color = [:orange :orange :orange :orange :blue :blue :blue :blue :blue :blue],bar_width=1,bins = xs)
savefig("all_neurons_conditions_histogram_no_label.png")


hist_plot = histogram(xs,[situation_1,situation_3_1,situation_2, situation_3_2,situation_3,situation_3_3,situation_4,situation_3_4,situation_5,situation_3_5],alpha=0.5,label=["1.5g Banana" "1.5g Banana" "0.3g Banana" "0.3g Banana" "0.2ml Juice" "0.2ml Juice" "0.5ml Juice" "0.5ml Juice" "0.9ml Juice" "0.9ml Juice"],layout=(5,2),xlims=(0,10),ylims=(0,6),xticks=xticks,size=(1000,800),xlabel="Spike Count in Window", ylabel="Number of Trials", color = [:orange :orange :orange :orange :blue :blue :blue :blue :blue :blue],bar_width=1,bins = xs)
savefig("all_neurons_conditions_histogram_juice_vs_standard.png")


age = 100*rand(1000,1)

histogram(clamp.(age, 19, 61), bins = 10; xlabel = "Age", leg = false,
                 xticks = ([17.5; 20:5:60; 62.5], ["<20"; 20:5:60; ">60"]), xrot = 45)
# there's no getting aorund this, you will have to make the global title in powerpoint or whatever you use for this

h = fit(Histogram, [situation_1,situation_2_1,situation_2, situation_2_2,situation_3,situation_2_3,situation_4,situation_2_4,situation_5,situation_2_5], layout=(5,2) nbins=10)
# fucking Rafal. He's completely ridiculous. Do I REALLY need to do this just so he can have some fucking BAR LABELS CENTERED -- taking up like an hour of time for FUCKING CENTERED BAR LABELS
#FUCK HIM I'm SO GLAD IM LEAVING!

using StatsBase
h = histogram(rand(100), nbins=10)
plot(h)
h.weights


xs2 = [1,2,3]
ys = [3,5,2]
histogram(xs2, ys,xticks=xs2)
b

## TODO -- do normalization by firing rate in period preceding cue and then difference in activity on the trial minus the baseline
# ensure that correlation between active neurons

# TODO normalize sequence of regressors by standard deviation -- subtract mean and divide by standard deviation for each sequence
# this is important because cannot reall interpret the beta weights express significance of factors and absolute value of regressor

## TODO how to see theory, try to assume -- encode what combination of two prediction errors each neuron encodes, if 30 neurons with each combination
# how to design a system that can accurately estimate mean quantity of juice vs banana
# i.e show for this particular combination and design system that learns this and simulate this system and produces same pattern of activity as in data

# for histogram figures, add label for 0, make sure labels are under the middles of the graph
# yellow by orange + thicker lines in the line graphs
# have moving window of 200 ms centered every 50 ms
# in histogram plot use same color scheme, for bananas make all bars orange
# asme thing for histogram to make sure labels are exactly in middle and in middle of bars
# ensure all bars have same width

#5 = 0.9ml juice
#1 = 1.5kg bananas
# 4 = 0.5ml juice
# 3 = 0.2ml juice
# 2  = 0.3kg bananas
#histogram!(xs,count_list_juice_neuron_food,nbins=10, alpha=0.5, label="Food Response")
# this is looking vaguely encouraign
spiketimes,  bits, situations, fp_times, stim_onsets, solenoid_1s, solenoid_2s, duration = parse_data_dict("Data/w065-0368" *".mat")
count_list_food_neuron_juice = counts_per_trial_window_by_situation(spiketimes, stim_onsets,situations, JUICE_SITUATIONS,100,400,false)
count_list_food_neuron_food = counts_per_trial_window_by_situation(spiketimes, stim_onsets, situations, BANANA_SITUATIONS, 100,400,false)
xs = collect(1:11)
histogram(xs,count_list_food_neuron_juice,nbins=10, alpha=0.5, label="Juice Response",bins=xs)
histogram!(xs ,count_list_food_neuron_food,nbins=10, alpha=0.5,label="Food Response")
xlabel!("Spike Count in Window")
ylabel!("Number of Trials")
title!("Food Neuron Response to Food/Juice")
savefig("food_neuron_response_histogram.jpg")


# okay, I think this is looking pretty nice. I'll then have to send this to Rafal and then be DONE with this
# then the only remaining thing to do is the fuckign PCN stuff. I'll have to leave the temporal PC stuff till tomorrow
# but can fiddle with some graphs and do some minor updates to get this working
# and it's nice to be able to do this
# somewhere I hasve the mean count graph. let's try to find this because if I can find this it could be quick also.
# these things don't take very long if you are focused but UGH staying focused is going to be tough
# maybe after you get back you can look at fertility test + verses so only thing left to do is STUPID temporal PC!
# although looks like tomorrow so let's ahve a look at it to figure it all out from scratch
# and have fun
bib
# TODO 1.) rerun shuffling over condition type and try to get that to work
# TODO more visualization, plot histogram for different volumes or water
# TODO do an armin style line plot distribution but keep colors close for banana/juice

# right let's do this next. First off I have to figure out how to do SHUFFLING
# let's FOCUS on this and get it done -- it IS highly achievable you just need to focus on executive functioning stuff
#figures, keep yaxis same, high quality,  -- don't export by jpeg, change by pdf or png
# move legend to side / keep only in one figure, remove title and set to mnist etc, remove caption in the thing
# cut down text if necessary -- UHN camera ready stuff


# show distribution across trials in a plot, i.e. over a trail, somek ind of histogram

# do dopamine neurons work over multiple days
# could be across multiple days, but monkeys can have fixed preferences

# so the result is highly significant given the bins
# okay great this is HIGHLY significant!! awesome! My result is super significant super exciting

# okay what else can be done to improve this analysis
# man this is taking SUCH a long time and you're being SO INEFFICIENT
# it's quite silly to be honest, at least you can make the histogram quite straightforward
# okay for 1000 going to be waiting for a LONG while... so there is that
# it's quite funny, I think these results will be okay
# and hopefully this histogram will be better at aroudn the correct time so we'll see
# but seems quite reasonable actually running like 1000 regressions
# I'm thinking this is pretty encouraging actually to see what's going on
# so it's actually significantly more and then I'll avoid the actual thing to see what's going on
# I have 10 / 19 significant reuslts which is actually pretty big
# wish I stored the histogram and it's actually quite important
# it's not totally convincing though as still having 7 be significant is quite big
!((500 > 100) && (500 < 400))
# this will take a while... so I'll need some basis level for significance
# i.e. less than 50 in 1000
# right so basically what I need to do is just define some kind of null hypothesis
# how can I do this... I can either do it by randomly selecting over window sizes
# all of size 300, but ignoring the current position, and then rerun the analysis
# or alternatively, I can shuffle assignment of situations?
# that's the alternative thing to do
# this is the last thing I need to do. let's randomly assign window sizes

# so it seems to be doing starts randomly, which isn't great
# I need to doible check that boolean logic as it's being annoying
# but would be interesting to look at it to see how much significant I'm actually getting
# very exciting though. and we might be under the 5% limit in which case it would seem to be significant
# it's an awful lot of experimental running though so ugh!
rand(Uniform(-8000,8000))
coeff_list
using NPZ
# I can always use the python scripts to properly plot the coefficients
# and see what is going on

# great so I've got these coefficient lists, what do I now do with them? it's hard to figure this out but doabel
# I can GET RESULTS SOON! so thati s important but UGH DOING SO BADLY FUCK

# right, let's figure out what I'm meant to do with these.  Also it really makes you realize that your salary is fucking good for london, better than even a lot of finance people even in their 30s
# i.e. you are at an elite technical level and from there you just need to level up over and over again
# and in fact that you can consider buying a STARTER HOME IN LONDON is fucking crazy
# that salary is doing okay actually it's in top 1% income just about hitting top 1% but that includes all people in fnance in like their 40s and 50s so it's doing okay overall
# to start to get some decent money to increase wealth by approx 6k/mo
# that is a damn slow increase to wealth however and no way to real wealth
# but you just need to FOCUS to get it sorted out
# okay got to focus to get this sorted. I think I actually fixed it and made it worse
# basically everything goes to juice, I need to run the LM and ACCEPT CONJECTURE
# and just UGH. I'm REALLY STRUGGLING TO STAY FOCUSED AND NOW AM TIRED UGH.
# I'm NOT ENJOYING THIS I've still got a lot of work left to do ugh and I've got to stay productive if at all possile
# to even be able to present anything to mark and everyone and actually MAKE THE SLIDES
# play FOCUS it's been like 1.5 hours and you have DONE NOTHING FUCK FUCK FUCK


bib
# right and baiscally for this neuron we want which has 43 datapoints, I guess that's not bad try to find least square single coefficients to solve them
# so there's just not thatm any counts. UGH. I really really hate this. but I enjoy coding and eventually I can get there too and try to figure it out effectively
# but it's justn ot working for me!
# what should I be regressioning -- you need to act FAST when you make a decision -- this is your other issue
# which is VERY BAD for investing and any kind of business is acting FAST and RESPONDING FAST TO EMAILS
# not doing so is obviously really stupid so who knows what's going on there an how well it is to do it etc
# and how to actually obtain correctness in all ways

# okay, let's try to respond fast to these things and actually do them reasonably. I REALLY hate it and I can just hand in my notice and not do any of the stupid Juan work
# to be honest although I shoudl do at least some of it to get it sorted out
# and have a whole bunch of deadlines and things to get done and then focus primarily on Yuhang and figuring out a few things relating to all my vairous projects to wind down
# let's try to sort this otu

# great well it looks like I already have the GLM code up there
# idk how exactly to determine significant but first I need to look at the results


# right so this isn't that useful what you want to do is figure out a linear regression across trials
# i.e. so basically what we want to do is to a.) code up the various situations
# and b.)
# basically we want to regress for each neuron, the number of spikes in the bin
# across each trial right?
# and for either amount or juice type
# and then see what the coefficients are and if they are significantly different from null
# i.e. we also need to get a subjective value estimate regressor we can read off from the paper
# and include all situations
#so the regression is just subjective value + Juce/banana
# it's essentially the same regression that I've already done, this is why it's so fucking STUPID that you procrastinated so badly



# right so this is showing same thing. Some neurons have like 80% of their spikes go to juice
# or at least above 70% which is quite big actually
# so yeah Idk what I'm seeing but like we need to normalize this a little
# but like yeah we just need to get some primarily responding to
mean_1s, mean_2s, stds_1s, stds_2s = neurons_spike_proportions(neuronlist, JUICE_SITUATIONS, BANANA_SITUATIONS,100,400)
mean_1s
scatter(mean_1s, mean_2s, legend=:none, xerr=stds_1s ./ sqrt(106), yerr = stds_2s ./sqrt(106))
xs = collect(0:0.1:5)
plot!(xs, xs)
xlabel!("Spikes Juice")
ylabel!("Spikes Banana")
title!("Juice spikes vs banana spikes by neurons")

L = length(neuronlist)


counts_juice, counts_banana = neuron_means_std_spikes_in_window_by_situation(neuronlist,[JUICE_SITUATIONS, BANANA_SITUATIONS],100,400,false, true)
counts_juice_means, counts_juice_stds = counts_juice
counts_banana_means, counts_banana_stds = counts_banana
counts_juice_means
#counts_banana_means = counts_banana_means ./ mean(counts_banana_means)
#counts_juice_means = counts_juice_means ./ mean(counts_juice_means)
# I'll have to write a function to split these apart. Super irritating and silly
# idk what is going on with that but I can do this all for rafal bullshit

# let's do a stupid function to do this
function new_filter_juice_banana(counts_juice, counts_banana,counts_juice_stds, counts_banana_stds)
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

for i in 1:length(counts_juice_means)
    println("hello")
end
counts_juice_stds
counts_banana_stds

# right yeah my exceptionally bad coding habits have really cocme home to make this increasingly unwieldy
# to try to actually figure out what is happening here
# I really think I need to move to the new method
# so respond to rafal that in fact this is the max
# and that this is what the other one looks like
# I majorly don't understand
# shit I also need to actually make the other figures which is super irritating
# dagnabbit rafal and he'll want to fiddle with this stuff
# it takes so long to debug and I don't know why
# super shit
# right yeah let's extract it all again and see if I can work with extracted verison
# as this is horrendous

#juice_more_juice = counts_juice_means[Int.(counts_juice_means .> counts_banana_means)]
scatter(counts_juice_means, counts_banana_means, legend=:none, xerr=counts_juice_stds ./ sqrt(106), yerr = counts_banana_stds ./sqrt(106),ylim=(-3.5,3.5),xlim=(-3.5,3.5),lc=:reds,mc=:reds, msc=:reds)
#scatter!([1,2,3],[2,2,2])
xs = collect(-3:0.1:3)
plot!(xs, xs,color=:black,linestyle=:dash,alpha=0.7)
xlabel!("Relative Spikes Juice")
ylabel!("Relative Spikes Banana")
title!("Juice spikes vs banana spikes by neurons")
savefig("figures/spikes_scatter_100_400_normalized.jpg")

# correlation coefficient
using Statistics
float_juice_counts = float.(counts_juice_means)
float_banana_counts = float.(counts_banana_means)
corr_coeff = Statistics.cor(float_juice_counts, float_banana_counts)
# R = 0.145


using HypothesisTests
p_val = pvalue(CorrelationTest(float_juice_counts, float_banana_counts))
# PVAL = 0.552
# scatter w/ diff colors
more_juice_juice, more_juice_banana,more_banana_juice,more_banana_banana,more_juice_juice_stds,more_juice_banana_stds,more_banana_juice_stds,more_banana_banana_stds = new_filter_juice_banana(counts_juice_means, counts_banana_means,counts_juice_stds, counts_banana_stds)
#scatter(more_juice_juice, more_juice_banana,xerr = more_juice_juice_stds ./ sqrt(106), yerr=more_juice_banana_stds ./sqrt(106),m=:circle,lc=:blue, mc=:blue, msc=:blue,alpha=1,ylim=(-3.5,3.5),xlim=(-3.5,3.5),label=false)
scatter(more_juice_juice, more_juice_banana,xerr = more_juice_juice_stds ./ sqrt(106), yerr=more_juice_banana_stds ./sqrt(106),m=:circle,lc=:blue, mc=:blue, msc=:blue,alpha=1,label=false)
scatter!(more_banana_juice, more_banana_banana,xerr = more_banana_juice_stds ./ sqrt(106), yerr=more_banana_banana_stds ./sqrt(106),m=:circle,lc=:orange, mc=:orange, msc=:orange,alpha=1,label=false)
scatter!(more_juice_juice, more_juice_banana,m=:circle,lc=:black, mc=:black, msc=:black,label=false)
scatter!(more_banana_juice, more_banana_banana,m=:black,lc=:black, mc=:black, msc=:black,label=false)
vline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
hline!([0], linestyle=:dash, alpha=0.6, color=:gray,label=false)
# right so what does this mean
xs = collect(-3:0.1:3)
plot!(xs, corr_coeff .* xs,label="Best fit line",color=:red,legend=:topleft)

#plot!(xs, xs,color=:black,linestyle=:dash,alpha=0.7,label=false)
xlabel!("Relative Spikes Juice")
ylabel!("Relative Spikes Banana")
title!("Juice spikes vs banana spikes by neurons")
savefig("figures/spikes_scatter_100_400_normalized_colored_correlation_no_best_fit.jpg")
# this is another issue, I'm just not interested in this so I do it very poorly
# and then it hurts; it hurts majorly so UGH
# quite irritating in its own way so we'll see waht's going on with that
# so yeah!


y1, y2 = rand(9), rand(9)
dy1, dy2 = 0.2rand(9), 0.2rand(9)
plot(y1, yerror=dy1, m=:square, lc=:reds, mc=:reds, msc=:reds)
plot!(y2, yerror=dy2, m=:square, lc=:blues, mc=:blues, msc=:blues)
scatter!(y1, m=:square, mc=:reds, label=false)
scatter!(y2, m=:square, mc=:blues, label=false)

b




# shit is there really much hope of me remmebering how to do linear regressions in julia for this?
# probably not but FUCK. I've been so busy with various things it's really really hurting and FUCK I just feel so sad and I'm so stressed
#okay THINK. It's not that hard. what do you have to do

# great so what am I meant to do here? It's really ough and I can't even remember?
# I'm meant to run a linear regression of SOMETHING
# and I think I did. ARGH. I FUCKING HATE THIS COMING BACK TO THIS AFTER SO LONG I JUST CANT BE BOTHERED
# and FUCK. JUANs stuff. I'm SO GLAD IM LEAVING AND NOT DOING THAT SOULSUCKING JOB
# I SUCK at doing things I hate and this is what I hate so thank god will be moving onto something different
# FUCK MY LIFE I'm just so stressed and not having fun at all
# FUCK
# so we'll see what happens here
# it will be quite interesting
# so right I mean there could be an effect here but it's hard to tell
# like there's definitely some neurons which are significantly higher for banana vs juice
# and vice versa but we don't really know and we are summing across trials
# lot's of neurons in the middle. Not perfect results
# so there's definitely a very strong effect of overall utility
# but what do we want to know about differential sensitivity -- we essentially have to subtract this component out right
# so I mean there's a couple of neurons on either side which are definitely far enough away
# even for the smaller ones there is some hope but i'm not totally sure what to make of this actually
# like what is going on here to figure it out and argh!
# rafal is really irritating to me and I've got him to submit to DM if he does email it me that will be huge
# DM is going w/ market so hopefully if it's equivalent comp at DM and VASTLY more stable you can really blossom and have so much fun there
# and do a startup in a few years that kind of thing basically
# and immediately start family planning so yeah by moving out you get all the negotiating leverage that you would want
# and it'll make a whole load of sense so that's good
# so yeah that was definitely worth going out in person to do things -- extremely worthwhile and very useful to chat about it
# personally, I feel that the reason that this is going on is that we haven't quite correctly normalized it
# I mean I'm not totally sure waht to do with these results
# like there are definitely -- realistically, I'm not totally sure I could tell with these data
# rafal says they are random but the error bars are good
# we just need some sensible way to normalize!

x = [1,2,3]
y = [4,5,6]
for (x_i,y_i) in zip(x,y)
    println("x: $x_i, y:$y_i")
end

# right, I've got to decide what to show rafal or not. It's tough and to be honest it's a bit irritating this
# but idk how to completely solve it
# it's very exciting re deepmind however.
# right so yeah I mean this isn't looking good overall
# I'm not totally sure what I should report but tbh I should just be honest about it



## next rafal idea TODO is to run a regression analysis using the subjective values extracted from the paper
# on the thing to see if there is any example and then do a shuffling thing to make sure that the number of significant regression coefficients is valid
# if we can believe this statistically then we would need to also simulate range of value of RB model essentially
# I'm pretty sure you can do anything with this and then into nature neuroscience
#
