import.pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
using CSV, DataFrames

dir = joinpath(@__DIR__, "..")
csv_file_path_final_demand_1 = joinpath(dir, "CAISO_Data", "Demand_Datafiles", "Final_Demand_1.csv")
csv_file_path_final_demand_2 = joinpath(dir, "CAISO_Data", "Demand_Datafiles", "Final_Demand_2.csv")
# Write the CSV file to the specified path

## Imort Net Demand for first year
Demand_1 = CSV.File(csv_file_path_final_demand_1) |> DataFrame
display(Demand_1)
# converting demand values from UTC to PST time (1am PST = 9 am UTC) changing column indexes
Demand_1 = Demand_1[!,[:9,:10,:11,:12,:13,:14,:15,:16,:17,:18,:19,:20,:21,:22,:23,:24,:1,:2,:3,:4,:5,:6,:7,:8]]
display(Demand_1)
Demand_1_array = Array{Float64, 2}(undef, 365, 24)
for i in 1:365
    for j in 1:24
        Demand_1_array[i,j] = Demand_1[i,j]
    end
end
display(Demand_1_array)

## Import Net Demand for second year
Demand_2 = CSV.File(csv_file_path_final_demand_2) |> DataFrame
# converting demand values from UTC to PST time (1am PST = 9 am UTC) changing column indexes
Demand_2 = Demand_2[!,[:9,:10,:11,:12,:13,:14,:15,:16,:17,:18,:19,:20,:21,:22,:23,:24,:1,:2,:3,:4,:5,:6,:7,:8]]
display(Demand_2)
Demand_2_array = Array{Float64, 2}(undef, 365, 24)
for i in 1:365
    for j in 1:24
        Demand_2_array[i,j] = Demand_2[i,j]
    end
end
display(Demand_2_array)


##clustering for net demand of first year
Final_Demand=hcat(Demand_1_array,Demand_2_array)

Demand_T = transpose(Final_Demand)
display(Demand_T)
Pkg.add("Clustering")

using Clustering

# cluster Net_Demand_T into 7 clusters using K-means
cluster = kmeans(Demand_T, 7; maxiter=1000, display=:iter)

@assert nclusters(cluster) == 7 # verify the number of clusters

a = assignments(cluster) # get the assignments of points to clusters
println(a)
counts(a)

## group of days for first year
daysGroup_1 = zeros(7,24)
# add all column vectors of each day from same group
for i in 1:365
        daysGroup_1[a[i],:] = daysGroup_1[a[i],:] + Demand_1_array[i,:]
end
display(daysGroup_1)
# devide each column vector by counts of each cluster
for i in 1:7
    daysGroup_1[i,:] = daysGroup_1[i,:] / counts(cluster)[i]
end

display(daysGroup_1)
##normalize demand data
d=findmax(daysGroup_1)

display(d[1])

## devide every value of the matrix by max value of demand data
Normalized_demand_1 = daysGroup_1 / d[1]
display(Normalized_demand_1)
Actual_Demand = 650*Normalized_demand_1
display(Actual_Demand)

d_bus = [0.1, 0 , 0.15, 0.2, 0.25, 0.05, 0, 0, 0.1, 0.15, 0]
demandData=Dict()

for i in 1:11
    demandData[string(i)]= Dict()
        for d in 1:7
            demandData[string(i)][string(d)] = Dict()
            for t in 1 : 24
                demandData[string(i)][string(d)][string(t)] = Dict("demand"=> d_bus[i]*Normalized_demand_1[d,t])
            end
        end
end

csv_file_path_solar_1 = joinpath(dir, "CAISO_Data", "Solar-Datafiles", "SOLAR_NP_SP_ZP_1_df.csv")
Solar_1 = CSV.File(csv_file_path_solar_1) |> DataFrame
display(Solar_1)
Solar_1 = Solar_1[!,[:9,:10,:11,:12,:13,:14,:15,:16,:17,:18,:19,:20,:21,:22,:23,:24,:1,:2,:3,:4,:5,:6,:7,:8]]
display(Solar_1)
Solar_1_array = Array{Float64, 2}(undef, 365, 24)
for i in 1:365
    for j in 1:24
        Solar_1_array[i,j] = Solar_1[i,j]
    end
end
display(Solar_1_array)

solarGroup = zeros(7,24)
# add all column vectors of each day from same group of solar data
for i in 1:365
        solarGroup[a[i],:] = solarGroup[a[i],:] + Solar_1_array[i,:]
end

# devide each column vector by counts of each cluster
for i in 1:7
    solarGroup[i,:] = solarGroup[i,:] / counts(cluster)[i]
end

display(solarGroup)
##max value from wind data
s=findmax(solarGroup)

display(s[1])

## devide every value of the matrix by max value of solar data
Normalized_solar= solarGroup / s[1]
print(Normalized_solar)
display(Normalized_solar[:,[12]])


## import Wind first year data

csv_file_path_wind_1 = joinpath(dir, "CAISO_Data", "Wind_Datafiles", "Wind_NP_SP_1.csv")
Wind_1 = CSV.File(csv_file_path_wind_1) |> DataFrame
display(Wind_1)
Wind_1 = Wind_1[!,[:9,:10,:11,:12,:13,:14,:15,:16,:17,:18,:19,:20,:21,:22,:23,:24,:1,:2,:3,:4,:5,:6,:7,:8]]
display(Wind_1)
Wind_1_array = Array{Float64, 2}(undef, 365, 24)
for i in 1:365
    for j in 1:24
        Wind_1_array[i,j] = Wind_1[i,j]
    end
end
display(Wind_1_array)

windGroup = zeros(7,24)
# add all column vectors of each day from same group of wind data
for i in 1:365
        windGroup[a[i],:] = windGroup[a[i],:] + Wind_1_array[i,:]
end

# devide each column vector by counts of each cluster
for i in 1:7
    windGroup[i,:] = windGroup[i,:] / counts(cluster)[i]
end

display(windGroup)

##max value from wind data to normalize the wind power values
w=findmax(windGroup)

display(w[1])

## devide every value of the matrix by max value of wind data
Normalized_wind = windGroup / w[1]
display(Normalized_wind)

##
windData=Dict()
    for d in 1:7
            windData[string(d)]= Dict()
        for t in 1 : 24
                windData[string(d)][string(t)] = Dict("wind_curve"=>Normalized_wind[d,t])
        end
    end

solarData=Dict()

    for d in 1:7
            solarData[string(d)]= Dict()
        for t in 1 : 24
                solarData[string(d)][string(t)] = Dict("solar_curve"=>Normalized_solar[d,t])
        end
    end


Pkg.add("PyPlot")
using PyPlot
x = 1:24
y1 = Normalized_solar[1,:]
y2 = Normalized_solar[2,:]
y3 = Normalized_solar[3,:]
y4 = Normalized_solar[4,:]
y5 = Normalized_solar[5,:]
y6 = Normalized_solar[6,:]
y7 = Normalized_solar[7,:]

fig, ax = subplots(figsize=(11, 5))
# Plot the data and set the legend labels
ax.plot(x, y1, label="DayType_1",lw=2, linestyle=":")
ax.plot(x, y2, label="DayType_2",lw=2,linestyle="--")
ax.plot(x, y3, label="DayType_3",lw=2, linestyle="-.")
ax.plot(x, y4, label="DayType_4",lw=2,linestyle=":")
ax.plot(x, y5, label="DayType_5",lw=2)
ax.plot(x, y6, label="DayType_6",lw=2)
ax.plot(x, y7, label="DayType_7",lw=2, linestyle=":")
# Create the legend with multiple columns
leg = ax.legend(ncol=4, loc="lower center", bbox_to_anchor=[0.5, -0.30], fontsize=10)
# Set the title and font size of the legend
suptitle("Solar Irradiance vs Time Of The Day", fontsize=16)
# Remove top and right spines
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)

# Set the position of the left spine to x=0
ax.spines["left"].set_position(("data", 0))

# Set the position of the bottom spine to y=0
ax.spines["bottom"].set_position(("data", 0))

# Set the linewidth of the x and y axes
ax.spines["bottom"].set_linewidth(1)
ax.spines["left"].set_linewidth(1)

ax.set_xlabel("Hours Of The Day", fontsize=14)
ax.set_ylabel("Solar Irradiance", fontsize=14)
ax.set_xticks(1:2:24)
# Show the plot
display(fig)
current_dir = @__DIR__

cd(joinpath(current_dir,"Input_Plots"))
savefig("Solar_Input_Plot.png", dpi=300, bbox_inches="tight")

## Input Plots for normalized wind output data of all day types

using PyPlot
x = 1:24
y1 = Normalized_wind[1,:]
y2 = Normalized_wind[2,:]
y3 = Normalized_wind[3,:]
y4 = Normalized_wind[4,:]
y5 = Normalized_wind[5,:]
y6 = Normalized_wind[6,:]
y7 = Normalized_wind[7,:]

fig, ax = subplots(figsize=(10, 5))
# Plot the data and set the legend labels
ax.plot(x, y1, label="DayType_1",lw=2, linestyle=":")
ax.plot(x, y2, label="DayType_2",lw=2,linestyle="--")
ax.plot(x, y3, label="DayType_3",lw=2, linestyle="-.")
ax.plot(x, y4, label="DayType_4",lw=2,linestyle=":")
ax.plot(x, y5, label="DayType_5",lw=2)
ax.plot(x, y6, label="DayType_6",lw=2)
ax.plot(x, y7, label="DayType_7",lw=2, linestyle=":")
# Create the legend with multiple columns
leg = ax.legend(ncol=4, loc="lower center", bbox_to_anchor=[0.5, -0.3], fontsize=10)
# Set the title and font size of the legend
suptitle("Normalized Wind Data vs Time Of The Day", fontsize=16)
# Remove top and right spines
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)

# Set the position of the left spine to x=0
ax.spines["left"].set_position(("data", 0))

# Set the position of the bottom spine to y=0
ax.spines["bottom"].set_position(("data", 0))

# Set the linewidth of the x and y axes
ax.spines["bottom"].set_linewidth(1)
ax.spines["left"].set_linewidth(1)
ax.set_xlabel("Hours Of The Day", fontsize=14)
ax.set_ylabel("Normalized Wind Data", fontsize=14)
ax.set_ylim(0,1)
ax.set_xticks(1:2:24)
# Show the plot
display(fig)

current_dir = @__DIR__

cd(joinpath(current_dir,"Input_Plots"))
savefig("Wind_Input_Plot.png", dpi=300, bbox_inches="tight")

##
using PyPlot
x = 1:24
y1 = Actual_Demand[1,:]
y2 = Actual_Demand[2,:]
y3 = Actual_Demand[3,:]
y4 = Actual_Demand[4,:]
y5 = Actual_Demand[5,:]
y6 = Actual_Demand[6,:]
y7 = Actual_Demand[7,:]
fig, ax = subplots(figsize=(7, 3))
# Plot the data and set the legend labels
ax.plot(x, y1, label="DayType_1",lw=2, linestyle=":")
ax.plot(x, y2, label="DayType_2",lw=2,linestyle="--")
ax.plot(x, y3, label="DayType_3",lw=2, linestyle="-.")
ax.plot(x, y4, label="DayType_4",lw=2,linestyle=":")
ax.plot(x, y5, label="DayType_5",lw=2)
ax.plot(x, y6, label="DayType_6",lw=2)
ax.plot(x, y7, label="DayType_7",lw=2, linestyle=":")
# Create the legend with multiple columns
leg = ax.legend(ncol=4, loc="lower center", bbox_to_anchor=[0.5, -0.5], fontsize=10)
# Set the title and font size of the legend
suptitle("Demand in KVA vs Time Of The Day", fontsize=13)
# Remove top and right spines
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)

# Set the position of the left spine to x=0
ax.spines["left"].set_position(("data", 0))

# Set the position of the bottom spine to y=0
ax.spines["bottom"].set_position(("data", 0))

# Set the linewidth of the x and y axes
ax.spines["bottom"].set_linewidth(1)
ax.spines["left"].set_linewidth(1)

ax.set_xlabel("Hours Of The Day", fontsize=11)
ax.set_ylabel("Actual Demand in KVA", fontsize=11)
ax.set_ylim(0,800)
ax.set_xticks(1:2:24)
# Show the plot
display(fig)
current_dir = @__DIR__
cd(joinpath(current_dir,"Input_Plots"))
savefig("Demand_Input_Plot.png", dpi=300, bbox_inches="tight")
