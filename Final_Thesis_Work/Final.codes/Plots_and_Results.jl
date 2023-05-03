include("C:\\Final_Julia_Work\\Final_Thesis_Work\\Final.codes\\Grid_Data.jl")
include("C:\\Final_Julia_Work\\Final_Thesis_Work\\Final.codes\\IncedenceMatrix.jl")
include("C:\\Final_Julia_Work\\Final_Thesis_Work\\Final.codes\\Clustering.jl")
include("C:\\Final_Julia_Work\\Final_Thesis_Work\\Final.codes\\GridModel.jl")

# storing value of solar power & wind power getting from solution
P_s=JuMP.value.(Ps)
print(1000*P_s[2,1,4,8])
P_w=JuMP.value.(Pw)

# storing values of demand served getting from the solution
Ped_s = JuMP.value.(P_e_d)

# storing values of charging and v2g power of EV's getting from the solution
Pch=JuMP.value.(P_ch)
Pv2g=JuMP.value.(P_V2G)

# storing value of solar power in an array for day 1, year 1 and time 1 to 24
Psd_1 = zeros(24)
for t in T
    Psd_1[t] = sum(1000*(P_s[s,1,4,t]) for s in S)
end
print(Psd_1)

print(Normalized_solar[2, :])
print(P_s[2, 1, 2, :])
println("")
println("")
println("")
println("")


for t in T
    print(t)
    print("\t     \t")
    print(round.((Ps_cap[2]*(solarData[string(4)][string(t)]["solar_curve"]));sigdigits=3))
    print("\t     \t")
    print(round.(P_s[2, 1, 4, t];sigdigits=3))
    print("\t     \t")
    print(round.((sum(1000*Pd_max[1]*demandData[string(i)][string(4)][string(t)]["demand"] for i in I));sigdigits=3))
    println("")
end

# storing value of required demand for day 1, year 1 and time 1 to 24
Ped_1 = zeros(24)
    for t in T
        Ped_1[t] = sum(1000*Pd_max[1]*demandData[string(i)][string(4)][string(t)]["demand"] for i in I)
    end
display(Ped_1)

# storing value of wind power for day 1, year 1 and time 1 to 24
Pwd_1 = zeros(24)
for t in T
    Pwd_1[t] = sum(1000*(P_w[w,1,4,t]) for w in W)
end
display(Pwd_1)

# storing value of demand served for day 1, year 1 and time 1 to 24
Pds_1 = zeros(24)
for t in T
    Pds_1[t] = sum(1000*(Ped_s[i,1,4,t]) for i in I)
end
display(Pds_1)

# storing value of charging power of EV for day 1, year 1 and time 1 to 24
Pchd_1 = zeros(24)
for t in T
    Pchd_1[t] = sum(1000*Pch[e,t,4,1] for e in E)
end
display(Pchd_1)

# storing value of v2g power for day 1, year 1 and time 1 to 24
Pv2g_1 = zeros(24)
for t in T
    Pv2g_1[t] = sum(1000*Pv2g[e,t,4,1] for e in E)
end
print(Pv2g_1)

######################################################################################################################################################
df = DataFrame(Demand_Required=Ped_1, Solar_served=Psd_1, Wind_Served= Pwd_1,V2G_served=Pv2g_1, Demand_Served =Pds_1)
print(df)
# plotting the results for each day type
using PyPlot
x = 1:24
y1 = Ped_1
y2 = Pwd_1
y3 = Psd_1
y4 = Pds_1
y5 = Pchd_1
y6 = Pv2g_1

fig, ax = subplots(figsize=(11, 5))
# Plot the data and set the legend labels
ax.plot(x, y1, label="Required_Demand",lw=2, linestyle=":")
ax.plot(x, y2, label="Wind_Power_Served",lw=2,linestyle="--")
ax.plot(x, y3, label="Solar_Power_Served",lw=2, linestyle="-.")
ax.plot(x, y4, label="Total_Demand_Served",lw=2)
ax.plot(x, y5, label="EV_charging_Power",lw=2)
ax.plot(x, y6, label="EV_V2G_Power",lw=2)

# Create the legend with multiple columns
leg = ax.legend(ncol=4, loc="lower center", bbox_to_anchor=[0.5, -0.30], fontsize=10)
# Set the title and font size of the legend
suptitle("Power in KVA vs Time Of The Day", fontsize=16)
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
ax.set_ylabel("Power in KVA", fontsize=14)
ax.set_xticks(1:2:24)
# Show the plot
display(fig)
cd("C:/Final_Julia_Work/Final_Thesis_Work/Final.codes/Output_Plots")
savefig("Output_Plot.png", dpi=300, bbox_inches="tight")
