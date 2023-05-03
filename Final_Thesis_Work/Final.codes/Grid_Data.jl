## This file defines all the Grid related data

## section: for defining sizes

n_years = 10             ## total number of years
n_solar = 7             ## total number of solar units
n_wind = 4              ## total number of wind units
n_converters = 2        ## total number of converters
n_time = 24             ## total number of hours in a day
n_days = 7              ## total number of types of similar days(viz. sunny days, hot days, windy days, cludy days, rainydays)
n_bus = 11              ## total number of busses
n_line = 12             ## total number of lines
n_EV = 15              ## total number of Electric vehicles
n_charging_stations =3 ## total number of charging stations
r =  0.03               ## interest rate on installation
M = 1000                ## high weight constant

## section: for defining collections


Y = collect(1:n_years)                           ##set of years
Yn = collect(2:n_years)                          ##set of years
S = collect(1:n_solar)                           ##set of solar generation units
W = collect(1:n_wind)                           ##Set of Wind generation units
C = collect(1:n_converters)                     ##Set of converters
T = collect(1:n_time)                           ##Set of time of the day
Tn = collect(2:n_time)                          ##Set of time for energy of EV
D = collect(1:n_days)                           ##Set of days of the year
I = collect(1:n_bus)                            ##Set of Buses
L = collect(1:n_line)                           ##Set of lines
E = collect(1:n_EV)                             ##Set of EVs
CH = collect(1:n_charging_stations)             ##Set of charging stations

## section: for defining limits for collections
Sbase = 1e6                ##MW
Vbase = 12e3                 ##kv
Zbase = 144               ##V*V/S
Demand1 = 650               ##KVA
Solar_Capacity = 80       ##KVA(for each)
Wind_Capacity = 100       ##KVA(for each)
Converter_Capacity = 250


Vmin = 0.9*ones(length(I))          ## minimum voltage of a bus
Vmax = 1.1*ones(length(I))          ## maximum voltage of a bus

Pw_cap =0.1*ones(length(W))        ## forecated active power of a wind

Ps_cap = 0.1*ones(length(S))        ## forecasted active power of a solar

Pc_min = zeros(length(C))         ## minimum active power of a converter
Pc_max =0.25*ones(length(C))         ## maximum active power of a converter
Qc_min = zeros(length(C))        ## minimum recative power of a converter
Qc_max =0.25*ones(length(C))            ## maximum reactive power of a converters
Pd_max = [0.65 0.665 0.68 0.695 0.71 0.726 0.74 0.75 0.77 0.79]  ##peak demand
Pf = 0.8
#Sl_min = zeros(length(L))            ## minimum apperant power of a line
K = 3e8    ## what is the unit?       ## penalty factor (10000$ per MWhr)
Pch_max =0.04* ones(length(E))       ## maximum charging power of an EV(100KW)

Pv2g_max =0.04*ones(length(E))        ## maximum V2G power delivered by an EV(150KW)
##Pdisch = 0.25*ones(length(E))       ## maximum discharging power of an EV(250KW)
EV_min = 0.01*ones(length(E))            ## minimum energy stored in an EV
EV_max = 0.1*ones(length(E))       ## maximum energy stored in an EV(100KW)
## section: for defining  system topology/configuration

solar_bus=[3 5 7 8 9 11 10]                                ## bus numbers to which solar units are connected (4 solar units)

wind_bus=[1 2 4 6]                                    ## bus numbers to which wind units are connected (2 wind units)

c_bus_ac = [3 5]                                    ## AC side of bus numbers to which converters are connected (2 converters)
c_bus_dc = [7 11]                                   ## DC side of bus numbers to which converters are connected (2 converters)
EV_bus=[8 8 8 9 9 9 10 10 10 10 8 8 9 9 10]         ## bus number to which charging station is connected

line_from_bus=[2 1 2 4 2 3 6 7 8 8 9 10]            ## bus numbers from which lines are connected (12 lines)
line_to_bus=[1 4 4 5 3 6 5 8 11 9 10 11]            ## bus numbers to which lines are connected (12 lines)

AC_DC=[0 0 0 0 0 0 1 1 1 1 1 1]                     ## if line AC then 0, for DC =1 (12 lines)


## section: for defining line parameters/ properties for power flow equations

SL_MAX=1000*2*[120,120,200,80,100,80,170,120,120,120,100,100]/Sbase         ## line max power
line_length=[100,120,160,180,120,80,140,270,350,200,300,200]*10         ## lenght of lines in meter (12 lines)
r=0.028/304.8                                                           ## ohm/1000ft *(1000ft/m) = ohm/meter (resi of line)
x=0.037/304.8                                                           ## ohm per meter
R=(r*line_length)/Zbase                                                 ## array of resistances of lines (12 lines)
X=(x*line_length)/Zbase                                                 ## array of reactances of lines  (12 lines)
g=zeros(length(L))                                                      ## arrays conductance of a line
b=zeros(length(L))                                                      ## array of suscaptance of a line

# this for loop calculates conductance and suscaptance
# based on line parameters defined above
for l in 1:length(L)
    g[l]=real(inv(R[l] + im * X[l]))
    b[l]=imag(inv(R[l] + im * X[l]))
end

# this is to define a conductance and suscaptance matrix for the lines defined in above topology
GG=zeros(length(I),length(I))
BB=zeros(length(I),length(I))
for l in L
    GG[line_from_bus[l],line_to_bus[l]] = -g[l]
    GG[line_to_bus[l],line_from_bus[l]] = -g[l]
    BB[line_from_bus[l],line_to_bus[l]] = -b[l]
    BB[line_to_bus[l],line_from_bus[l]] = -b[l]
end
for i in I
    GG[i, i] = -sum(GG[i, :])
    BB[i, i] = -sum(BB[i, :])
end

## section: for defining constants used in objective function
cost_solar = 200000  # cost of 150kva solar panel
cost_wind = 580000  # cost of 200kva wind unit
cost_converter = 57000 # cost of 250kva converter unit
cost_line = 700000 # cost of line installation $7000 for per meter
cost_ChargingStation = 30000 # cost of 100kva charging station unit
