# Running the code - example:
# D:\julia-1.6.1\bin\julia.exe --color=yes --project=C:\Users\cesin\.julia\environments\v1.6 d:\GitRoot\Julia_InfectionSimulation\infection_simulation.jl --initial_vacc_rate 0.25 -d 500

# Running from the REPL:
# include("infection_simulation.jl")
# df = create_db_population()

# Import libraries

using Dates          # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("Dates")'
using ArgParse       # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("ArgParse")'
using Random         # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("Random")'
using Statistics     # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("Statistics")'
using DataFrames     # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("DataFrames")'
using Distributions  # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("Distributions")'
using StatsBase      # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("StatsBase")'
using CSV            # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("CSV")'
using Printf         # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("Printf")'
using Plots          # You may need to run in your command prompt: julia -e 'import Pkg; Pkg.add("Plots")'

function parse_commandline()
    # Retrieve command line arguments - see https://carlobaldassi.github.io/ArgParse.jl/latest/arg_table.html
    settings = ArgParseSettings()
    @add_arg_table settings begin
        "--num_trials", "-t"
            arg_type = Int64
            default = 1000
        "--max_iteration_days", "-i"
            arg_type = Int64
            default = 50
        "--population", "-p"
            arg_type = Int64
            default = 900
        "--max_days_sick", "-m"
            arg_type = Int64
            default = 12
        "--gender_prob", "-g"
            arg_type = Float64
            default = 0.52
        "--age_avg", "-a"
            arg_type = Int64
            default = 46
        "--age_std", "-s"
            arg_type = Int64
            default = 13
        "--comorbid_pois", "-c"
            arg_type = Float64
            default = 0.70
        "--asym_prob", "-y"
            arg_type = Float64
            default = 0.02
        "--self_isolate_days", "-l"
            arg_type = Int64
            default = 8
        "--antivax_prob", "-x"
            arg_type = Float64
            default = 0.09
        "--vax_effective_rate", "-e"
            arg_type = Float64
            default = 0.80
        "--initsick_prob", "-k"
            arg_type = Float64
            default = 0.1
        "--dayssic_default", "-d"
            arg_type = Int64
            default = 5
        "--initial_vacc_rate", "-v"
            arg_type = Float64
            default = 0.5
        #############################################################################################
        # The number that serves as the minimum for a random number (to max = 1.0)
        # above which people die if their probability of death > the random number
        "--death_threshold", "-o"
            arg_type = Float64
            default = 0.30
        #############################################################################################
        # Use 1.0 to make vaccination guarantted to save the sick from dying
        # or a number less than 1.0, just to decrease chances of dying
        "--vacc_saving_throw", "-w"
            arg_type = Float64
            default = 1.00
        #############################################################################################
        # R-naught value (each sick person can infect R0 healthy people per interaction i.e. per day
        # - not the traditional R0 which is not per interaction, but rather per infection)
        # This number has been carefully tuned and probably shouldn't be modified or if so, by little            
        "--r_naught", "-z"
            arg_type = Float64
            default = 0.1
        #############################################################################################
    end

    return parse_args(settings)
end

function create_db_population(
    population         = 900
  , day                = 0
  , gender_prob        = 0.52
  , age_avg            = 46
  , age_std            = 13
  , comorbid_pois      = 0.70
  , asym_prob          = 0.02
  , antivax_prob       = 0.09
  , initial_vacc_rate  = 0.50
  , vax_effective_rate = 0.80
  , initsick_prob      = 0.1
  , dayssic_default    = 5
)

    # https://dataframes.juliadata.org/stable/

    # Create an empty dataframe
    df = DataFrame()

    # https://discourse.julialang.org/t/add-columns-to-dataframe/49305

    # Generate an ID row having a sequence of consecutive numbers from 1 to population (# of rows)
    # https://discourse.julialang.org/t/how-to-generate-a-sequence-of-numbers-in-julia/19864
    df[!, "id"] = 1:1:population

    # Add the day number for each day in the trial
    df[!, "day"] .= day

    # Use Bernoulli distribution with probability gender_prob to generate true (Male) or false (Female), and repeat for each row in population
    # https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.Bernoulli
    df[!, "gender"] = rand(Distributions.Bernoulli(gender_prob), population)

    # Use Normal distribution with mean = age_avg, standard deviation = age_std and repeat for each row in population
    # https://juliastats.org/Distributions.jl/stable/starting/
    df[!, "Age"] = floor.(Int, round.(rand(Distributions.Normal(age_avg, age_std), population)))

    # Use Poisson distribution to generate the number of comorbidities for each subject in population
    # https://juliastats.org/Distributions.jl/v0.14/univariate.html#Distributions.Poisson
    df[!, "Comorbidity"] = rand(Distributions.Poisson(comorbid_pois), population)

    # Use Bernoulli distribution with probability asym_prob to generate true (Asymptomatic) or false (Symptomatic), and repeat for each row in population
    df[!, "Asymptomatic"] = rand(Distributions.Bernoulli(asym_prob), population)

    # Use Bernoulli distribution with probability anti_prob to generate true (Anti-vax) or false (Pro-vax), and repeat for each row in population
    df[!, "Antivax"] = rand(Distributions.Bernoulli(antivax_prob), population)

    # Next, determine who will get vaccinated
    # It will be a random sample of size of those who are pro-vax * vaccination rate
    # https://juliastats.org/StatsBase.jl/stable/sampling/
    vaccinated_ids = sample(df[(df.Antivax .== false), :id], floor.(Int, round(initial_vacc_rate * count(df.Antivax .== false))), replace=false, ordered=true)
    # Assume by default that the people are not vaccinated
    df[!, "Vaccinated"] .= false
    # Now update those who did gete vaccinated to true (some portion of those who are pro-vax)
    df[[x in vaccinated_ids for x in df[!, :id]], :Vaccinated] .= true

    # To determine whether the vaccine was effective, apply the probability of effectivenes of the vaccine to those who got vaccinated
    # Now draw the sample (without replacement) from a sample size based on vaccine effectiveness * those who got vaccinated
    #@printf("length(vaccinated_ids) = %d\n", length(vaccinated_ids))
    #@printf("floor.(Int, round(vax_effective_rate * length(vaccinated_ids))) = %d\n", floor.(Int, round(vax_effective_rate * length(vaccinated_ids))))
    vax_effective_ids = sample(vaccinated_ids, floor.(Int, round(vax_effective_rate * length(vaccinated_ids))), replace=false, ordered=true)
    # Assume by default that the vaccine is not effective (this will be true for anti-vax crowd AND pro-vax, for now)
    df[!, "Veffective"] .= false
    # Now set all pro-vax who were sampled to have an effective vaccination to true
    df[[x in vax_effective_ids for x in df[!, :id]], :Veffective] .= true

    # Next, anyone who has not been effectively vaccinated is subject to being infected at the start of the simulation as given by the parameter initsick_prob
    vax_ineffective_ids = sample(df[(df.Veffective .== false), :id], floor.(Int, round(initsick_prob * count(df.Veffective .== false))), replace=false, ordered=true)
    
    # Set the default of the number of days sick to 0
    df[!, "DaysSick"] .= 0
    # For those that started the simulation infected, have them start with the default number of days infected
    df[[x in vax_ineffective_ids for x in df[!, :id]], :DaysSick] .= dayssic_default

    # Set the default health status to 0 (0=never infected, 1=currently infected, 2=recovered, 3=deceased)
    df[!, "HealthStatus"] .= 0
    # For those that started the simulation infected, set their health status to 1
    df[[x in vax_ineffective_ids for x in df[!, :id]], :HealthStatus] .= 1

    # We can view the results of our dataframe creation using:
    # show(df, allrows=true, allcols=true, summary=true)

    return df

end

function create_db_results()

    df = DataFrame(
          TrialID               = Int[]
        , Day                   = Int[]
        , CountVaccinated       = Int[]
        , CountEffectiveVaccine = Int[]
        , CountNeverInfected    = Int[]
        , CountSick             = Int[]
        , CountRecovered        = Int[]
        , CountDeceased         = Int[]
        , CountCumSick          = Int[]
        , CountCumRecovered     = Int[]
        , CountCumDeceased      = Int[]
    )

    return df

end

function configure_text_overlay(
    ; text_loc_x              =   0.5
    , text_loc_y              =   0.5
    , text                    =   "CESAR\nMUGNATTO"
    , alignment               =   :center
    , fontfamily              =   "Verdana"
    , fontscale               =   0.1
    #, color                   =   RGBA(0, 1, 0, 64.0/255.0) # Corresponds to #00FF0040
    , color                   =   RGBA(128.0/255.0, 0, 1, 96.0/255.0) # Corresponds to #8000FF60
    , rotation                =   45
)
    text_overlay = Dict(
          "text_loc_x"          => text_loc_x
        , "text_loc_y"          => text_loc_y
        , "text"                => text
        , "alignment"           => alignment
        , "fontfamily"          => fontfamily
        , "fontscale"           => fontscale
        , "color"               => color
        , "rotation"            => rotation
    )

    return text_overlay

end

function make_plots(
      df                   :: DataFrame
    , text_overlay         :: Union{ Dict, Nothing }
    , out_file_name        :: String
)

    daily_stats_line_plot(df, text_overlay, out_file_name)
    cumul_stats_line_plot(df, text_overlay, out_file_name)

end

function daily_stats_line_plot(
      df                   :: DataFrame
    , text_overlay         :: Union{ Dict, Nothing }
    , out_file_name        :: String
)

    # https://dataframes.juliadata.org/stable/man/working_with_dataframes/#Taking-a-Subset
    df_plot = df[!, [:Day, :CountSick, :CountRecovered, :CountDeceased]]
    # Alternatively:
    #df_plot = select(df, ["Day", "CountSick", "CountRecovered", "CountDeceased"], copycols=false)
    df_line_plot = combine(groupby(df_plot, ["Day"]), :CountSick => mean, :CountRecovered => mean, :CountDeceased => mean; renamecols=false)

    x_lim_min = 0
    x_lim_max = findmax(df_line_plot.Day)[1]
    y_lim_min = 0
    y_lim_max = findmax(df_line_plot.CountSick)[1]
    x_plot_size = 600
    y_plot_size = 600
    plot_title = "Mean daily rate of infected, recovered, deceased"
    x_label = "Day"
    y_label = "Count"
    labels = ["Infected" "Recovered" "Deceased"]
    colors = [:red :blue :black]
    plot_lin = plot(
        (df_line_plot.Day, [df_line_plot.CountSick, df_line_plot.CountRecovered, df_line_plot.CountDeceased])
        , size          = (x_plot_size, y_plot_size)
        , title         = plot_title
        , xlabel        = x_label
        , ylabel        = y_label
        , linewidth     = 3
        , linecolor     = colors
        , label         = labels
        , legend        = :topright # :right, :left, :top, :bottom, :inside, :best, :legend, :topright, :topleft, :bottomleft, :bottomright
    )
    if isa(text_overlay, Dict)
        text_loc_x = Int16(round((x_lim_max - x_lim_min) * text_overlay["text_loc_x"])) + x_lim_min
        text_loc_y = Int16(round((y_lim_max - y_lim_min) * (1- text_overlay["text_loc_y"]))) + y_lim_min
        fontsize = Int16(round(min(x_plot_size, y_plot_size) * text_overlay["fontscale"]))
        #println("($(text_loc_x), $(text_loc_y))")
        annotate!([
            (
                  text_loc_x
                , text_loc_y
                , Plots.text(
                    text_overlay["text"]
                    , text_overlay["color"]
                    , text_overlay["alignment"]
                    , fontsize
                    , text_overlay["fontfamily"]
                    , rotation=text_overlay["rotation"]
                )
            )
        ])
    end

    savefig("$(out_file_name)_daily.png")

end

function cumul_stats_line_plot(
      df                   :: DataFrame
    , text_overlay         :: Union{ Dict, Nothing }
    , out_file_name        :: String
)

    df_plot = df[!, [:Day, :CountCumSick, :CountCumRecovered, :CountCumDeceased]]
    # Alternatively:
    #df_plot = select(df, ["Day", "CountCumSick", "CountCumRecovered", "CountCumDeceased"], copycols=false)
    df_line_plot = combine(groupby(df_plot, ["Day"]), :CountCumSick => mean, :CountCumRecovered => mean, :CountCumDeceased => mean; renamecols=false)

    x_lim_min = 0
    x_lim_max = findmax(df_line_plot.Day)[1]
    y_lim_min = 0
    y_lim_max = findmax(df_line_plot.CountCumSick)[1]
    x_plot_size = 600
    y_plot_size = 600
    plot_title = "Mean cumulative infected, recovered, deceased"
    x_label = "Day"
    y_label = "Count"
    labels = ["Infected" "Recovered" "Deceased"]
    colors = [:red :blue :black]
    plot_lin = plot(
        (df_line_plot.Day, [df_line_plot.CountCumSick, df_line_plot.CountCumRecovered, df_line_plot.CountCumDeceased])
        , size          = (x_plot_size, y_plot_size)
        , title         = plot_title
        , xlabel        = x_label
        , ylabel        = y_label
        , linewidth     = 3
        , linecolor     = colors
        , label         = labels
        , legend        = :bottomright # :right, :left, :top, :bottom, :inside, :best, :legend, :topright, :topleft, :bottomleft, :bottomright
    )
    if isa(text_overlay, Dict)
        text_loc_x = Int16(round((x_lim_max - x_lim_min) * text_overlay["text_loc_x"])) + x_lim_min
        text_loc_y = Int16(round((y_lim_max - y_lim_min) * (1- text_overlay["text_loc_y"]))) + y_lim_min
        fontsize = Int16(round(min(x_plot_size, y_plot_size) * text_overlay["fontscale"]))
        annotate!([
            (
                text_loc_x
                , text_loc_y
                , Plots.text(
                    text_overlay["text"]
                    , text_overlay["color"]
                    , text_overlay["alignment"]
                    , fontsize
                    , text_overlay["fontfamily"]
                    , rotation=text_overlay["rotation"]
                )
            )
        ])
    end

    savefig("$(out_file_name)_cumul.png")

end

function main()

    dt_start = Dates.now()

    # Read the arguments passed in and store their values
    parsed_args = parse_commandline()

    # println(parsed_args)

    num_trials = parsed_args["num_trials"]
    max_iteration_days = parsed_args["max_iteration_days"]
    population = parsed_args["population"]
    max_days_sick = parsed_args["max_days_sick"]
    gender_prob = parsed_args["gender_prob"]
    age_avg = parsed_args["age_avg"]
    age_std = parsed_args["age_std"]
    comorbid_pois = parsed_args["comorbid_pois"]
    asym_prob = parsed_args["asym_prob"]
    self_isolate_days = parsed_args["self_isolate_days"]
    antivax_prob = parsed_args["antivax_prob"]
    vax_effective_rate = parsed_args["vax_effective_rate"]
    initsick_prob = parsed_args["initsick_prob"]
    dayssic_default = parsed_args["dayssic_default"]
    initial_vacc_rate = parsed_args["initial_vacc_rate"]
    death_threshold = parsed_args["death_threshold"]
    vacc_saving_throw = parsed_args["vacc_saving_throw"]
    r_naught = parsed_args["r_naught"]

    @assert (num_trials > 0) && (num_trials <= 2000) "Parameter 'num_trials' must be between 1 and 2000"
    @assert (max_iteration_days > 0) && (max_iteration_days <= 200) "Parameter 'max_iteration_days' must be between 1 and 200"
    @assert (population > 0) && (population <= 1000) "Parameter 'population' must be between 1 and 1000"
    @assert (max_days_sick > 0) && (max_days_sick <= 20) "Parameter 'max_days_sick' must be between 1 and 20"
    @assert (gender_prob >= 0.0) && (gender_prob <= 1.0) "Parameter 'gender_prob' must be between 0.0 and 1.0"
    @assert (age_avg > 20) && (age_avg < 100) "Parameter 'age_avg' must be between 21 and 99"
    @assert (age_std > 4) && (age_std < 20) "Parameter 'age_std' must be between 5 and 19"
    @assert (comorbid_pois > 0.0) && (comorbid_pois <= 5.0) "Parameter 'comorbid_pois' must be greater than 0.0 and less than or = 5.0"
    @assert (asym_prob >= 0.0) && (asym_prob <= 1.0) "Parameter 'asym_prob' must be between 0.0 and 1.0"
    @assert (self_isolate_days >= 0) && (self_isolate_days <= max_days_sick) "Parameter 'self_isolate_days' must be between 0.0 and max_days_sick."
    @assert (antivax_prob >= 0.0) && (antivax_prob <= 1.0) "Parameter 'antivax_prob' must be between 0.0 and 1.0"
    @assert (vax_effective_rate >= 0.0) && (vax_effective_rate <= 1.0) "Parameter 'vax_effective_rate' must be between 0.0 and 1.0"
    @assert (initsick_prob >= 0.0) && (initsick_prob <= 1.0) "Parameter 'initsick_prob' must be between 0.0 and 1.0"
    @assert (dayssic_default > 0) && (dayssic_default <= 8) "Parameter 'dayssic_default' must be between 1 and 8"
    @assert (initial_vacc_rate >= 0.0) && (initial_vacc_rate <= 1.0) "Parameter 'initial_vacc_rate' must be between 0.0 and 1.0"
    @assert (death_threshold >= 0.0) && (death_threshold <= 1.0) "Parameter 'death_threshold' must be between 0.0 and 1.0"
    @assert (vacc_saving_throw >= 0.0) && (vacc_saving_throw <= 1.0) "Parameter 'vacc_saving_throw' must be between 0.0 and 1.0"
    @assert (r_naught >= 0.0) && (r_naught <= 1.0) "Parameter 'r_naught' must be between 0.0 and 1.0"

    # For debugging (1)
    #return parsed_args

    # For debugging (2)
    #return create_db_population(
    #    population
    #  , gender_prob
    #  , age_avg
    #  , age_std
    #  , comorbid_pois
    #  , asym_prob
    #  , antivax_prob
    #  , initial_vacc_rate
    #  , vax_effective_rate
    #  , initsick_prob
    #  , dayssic_default
    #)

    gender_dict = Dict(
          1 => "M"
        , 0 => "F"
    )
    health_dict = Dict(
          0 => "Never infected"
        , 1 => "Sick"
        , 2 => "Recovered"
        , 3 => "Deceased"
    )

    # Ensure our output directory exists
    mkpath("output")

    # Create a new results database
    df_trials = create_db_results()

    # Create a matrix of zeros sized as the number of trials x number of days
    trial_results = zeros(Int64, num_trials, max_iteration_days)

    # Iterate through each trial
    for n = 1:num_trials

        # For each trial, we create a new population database
        # Create a new randomized population for this trial
        df_sim = create_db_population(
              population
            , 0
            , gender_prob
            , age_avg
            , age_std
            , comorbid_pois
            , asym_prob
            , antivax_prob
            , initial_vacc_rate
            , vax_effective_rate
            , initsick_prob
            , dayssic_default
        )

        count_vaccinated = length(df_sim[(df_sim.Vaccinated .== true), :id])
        count_effective_vaccine = length(df_sim[(df_sim.Veffective .== true), :id])

        # Loop over the maximum number of days we are testing
        for i = 1:max_iteration_days

            df_sim.day .= i

            # For each day in the trial, we need to:
            # 1. Determine which population members are infectious
            #   a. Must be in status "sick"
            #   b. must not have quarantined (quarantine is affected by asymptomatic status)
            # 2. Determine which population members can become sick
            #   a. Must be currently in a "Never been sick" status
            #   b. Must be unvaccinated (AntiVax) or vaccinated, but with a failed vaccination
            # 3. Apply the R0 (R-naught) value selected for the infection
            #    to see how many of the healthy are infected by the sick

            # Get the ids of people who are sick (HealthStatus==1)
            # AND either are not self-isolating, or are asymptomatic
            ids_sick = df_sim[
                (df_sim.HealthStatus .== 1) .& ((df_sim.DaysSick .< self_isolate_days) .| (df_sim.Asymptomatic .== true))
                , :id
            ]

            # Get the ids of people who are healthy (HealthStatus==0)
            # AND either have effective vaccination, or are AntiVax
            ids_healthy = df_sim[
                (df_sim.HealthStatus .== 0) .& (df_sim.Veffective .== false)
                , :id
            ]

            # Multiply the number of currently infectious people by R0 to determine
            # how many of the currently healthy people will get infected
            # We choose the minimum of the number of healthy people remaining
            # vs the number of sick (who have not self-isolated) * R0
            # This is because with replace=False, we cannot choose
            # more than the number of healthy people remaining
            # Use the value of r_naught as the upper limit to choose a
            # random uniform number of new infections for the current day
            day_infection_rate = rand(Distributions.Uniform(0.0, r_naught))
            ids_newly_infected = sample(ids_healthy, min(length(ids_healthy), floor.(Int, round(length(ids_sick) * day_infection_rate))), replace=false, ordered=true)

            # Update the population database
            # 1. Add the newly-infected
            df_sim[[x in ids_newly_infected for x in df_sim[!, :id]], :HealthStatus] .= 1

            # 2. Increment the sick days (up to the maximum)
            # for everyone currently infected, including the newly-infected
            ids_updated_sick = df_sim[(df_sim.HealthStatus .== 1), :id]
            # Add 1 to each member in the array for the current days sick, but do not exceed the max_days_sick
            df_sim[[x in ids_updated_sick for x in df_sim[!, :id]], :DaysSick] .= min.(
                df_sim[[x in ids_updated_sick for x in df_sim[!, :id]], :DaysSick] .+ 1, max_days_sick
            )

            # 3. Determine recovery or death status if maximum days have been reached
            ids_sick_ended = df_sim[df_sim.DaysSick .>= max_days_sick, :id]
            # and reset the number of days sick to 0
            df_sim[[x in ids_sick_ended for x in df_sim[!, :id]], :DaysSick] .= 0

            # 4. Determine who might have died from the virus
            # - use age, comorbidity as factors that contribute to likelihood of death
            # This seems a bit arbitrary but via experimentation it works
            # Take the age and divide by 200, then use this quotient as the exponent for Euller's number
            # The formula can be written as below by typing \euler (then pressing tab) - don't forget the dot operator for vector operations
            # prob_death_by_age = ℯ^(df_sim[[x in ids_sick_ended for x in df_sim[!, :id]], :Age] ./ 200.0 .- 1.0)
            # Alternatively, you can use:
            prob_death_by_age = exp.(df_sim[[x in ids_sick_ended for x in df_sim[!, :id]], :Age] ./ 200.0 .- 1.0)

            prob_comorbid_mult = exp.(df_sim[[x in ids_sick_ended for x in df_sim[!, :id]], :Comorbidity] ./ 20.0)

            # The following will assign a multiplier of 1
            # (no change to the above multipliers) for people who had not taken the vaccine
            # or some number less than 1 (reducing the chances of death)
            # for those who took the vaccine, but where it was not effective
            prob_death_vaccine_mult = (
                1 .- df_sim[[x in ids_sick_ended for x in df_sim[!, :id]], :Vaccinated] .* vacc_saving_throw
            )

            # Calculate the final probability of death based on age and comorbidities as factors
            # that contribute to likelihood of death and on the probability of the vaccine to prevent
            # deaths even if ineffective for preventing getting sick in the first place
            prob_death = round.(prob_death_by_age .* prob_comorbid_mult .* prob_death_vaccine_mult, digits=2)

            # Randomly choose a threshold for each person
            # to see if they recover or die based on their final probability
            rand_death_threshold = rand(Distributions.Uniform(death_threshold, 1.0), length(prob_death))
            # Get a vector of true/false values where the probability of death exceeds the death threshold for each patient
            actually_dead = prob_death .>= rand_death_threshold
            # Based on these true/false values, get the id values for those that are true
            idx_actually_dead = ids_sick_ended[actually_dead]
            # Those that are false are for people who did not die i.e. they recovered
            idx_recovered = ids_sick_ended[actually_dead .== false]

            df_sim[[x in idx_actually_dead for x in df_sim[!, :id]], :HealthStatus] .= 3 # Deceased
            df_sim[[x in idx_recovered for x in df_sim[!, :id]], :HealthStatus] .= 2 # Recovered

            #if n < 2
            #    outfile = @sprintf("output/trial(%04d)_%s.csv", n, Dates.format(dt_start, "yyyymmdd_HHMMSS"))
            #    CSV.write(outfile, df_sim, writeheader=(i==1), append=(i>1))
            #end

            # At the end of the day, append a new row to the results DataFrame
            # https://stackoverflow.com/questions/63100620/any-better-equivalent-for-pandas-value-counts-in-julia-dataframes
            dict_vals = countmap(df_sim.HealthStatus)

            # https://docs.julialang.org/en/v1/base/collections/#Base.get!
            CountNeverInfected = get!(dict_vals, 0, 0)  # Healthy
            CountSick          = get!(dict_vals, 1, 0)  # Sick
            CountRecovered     = get!(dict_vals, 2, 0)  # Recovered
            CountDeceased      = get!(dict_vals, 3, 0)  # Dead

            lst_row = [
                  n
                , i
                , count_vaccinated
                , count_effective_vaccine
                , CountNeverInfected                      # Healthy
                , CountSick                               # Sick
                , length(idx_recovered)                   # Recovered
                , length(idx_actually_dead)               # Dead
                , population - CountNeverInfected         # Cumulative Sick
                , CountRecovered # (already cumulative)   # Cumulative Sick Recovered
                , CountDeceased # (already cumulative)    # Cumulative Sick Dead
            ]

            push!(df_trials, lst_row)

            # Also, check if we can interrupt the looping now.
            # If the number of sick people still circulating has fallen to 0
            # then we no longer have to loop - every subsequent day will have the same results as the current day

            if length(df_sim[df_sim.DaysSick .== 0, :id]) == population #nrow(df_sim)
                break
            end

        end

    end

    dt_end = Dates.now()

    # https://dataframes.juliadata.org/stable/man/importing_and_exporting/

    #list_args = [[k, v] for (k,v) in parsed_args]
    #push!(list_args, ["start_datetime", Dates.format(dt_start, "yyyy-mm-dd HH:MM:SS")])
    #push!(list_args, ["end_datetime", Dates.format(dt_end, "yyyy-mm-dd HH:MM:SS")])

    outfileCSV = @sprintf("output/trials_%s.csv", Dates.format(dt_start, "yyyymmdd_HHMMSS"))
    outfilePNG = @sprintf("output/trials_%s", Dates.format(dt_start, "yyyymmdd_HHMMSS"))

    dict_params = Dict([(k,string(v)) for (k,v) in parsed_args])
    CSV.write(outfileCSV, sort(collect(dict_params), by=x->x[1]), delim=" = ", writeheader=false)
    
    dict_dates = Dict("start_datetime" => Dates.format(dt_start, "yyyy-mm-dd HH:MM:SS"), "end_datetime" => Dates.format(dt_end, "yyyy-mm-dd HH:MM:SS"))
    CSV.write(outfileCSV, sort(collect(dict_dates), by=x->x[2]), delim=" = ", writeheader=false, append=true)

    CSV.write(outfileCSV, df_trials, writeheader=true, append=true)

    text_overlay = configure_text_overlay()
    text_overlay = configure_text_overlay(
        ; text_loc_x              =   0.05
        , text_loc_y              =   0.5
        , text                    =   join([@sprintf("%s = %s", k, v) for (k,v) in sort(collect(dict_params), by=x->x[1])], "\n")
        , alignment               =   :left
        , fontfamily              =   "Verdana"
        , fontscale               =   0.02
        , color                   =   RGBA(128.0/255.0, 0, 1, 96.0/255.0) # Corresponds to #8000FF60
        , rotation                =   0
    )

    #make_plots(df_trials, nothing, outfilePNG)
    make_plots(df_trials, text_overlay, outfilePNG)

end

main()