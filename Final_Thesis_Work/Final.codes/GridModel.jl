cd("C:\\My_Julia_Work_WB_1\\My_Thesis_Work")

Pkg.resolve()
pwd()
Pkg.activate(".")
Pkg.instantiate()
Pkg.status()

import Pkg
Pkg.add("JuMP")
Pkg.add("Mosek")
Pkg.add("Gurobi")
using JuMP, Gurobi

## including grid data file which defines all the system and network parameter
#include("C:\\Final_Julia_Work\\Final_Thesis_Work\\Final.codes\\Grid_Data.jl")
#include("C:\\Final_Julia_Work\\Final_Thesis_Work\\Final.codes\\IncedenceMatrix.jl")
#include("C:\\Final_Julia_Work\\Final_Thesis_Work\\Final.codes\\Clustering.jl")
## defining model for optimiziation
GridModel = Model(Gurobi.Optimizer)

## defining model variables
@variable(GridModel, Is[s in S,y in Y] ,Bin )                                   ## binary variable solar installation
@variable(GridModel, Iw[w in W,y in Y] ,Bin )                                   ## binary variable wind installation
@variable(GridModel, Ic[c in C,y in Y] ,Bin )                                   ## binary variable converter installation
@variable(GridModel, Il[l in L,y in Y] ,Bin )                                   ## binary variable line installation


@variable(GridModel, V[i in I, y in Y, d in D, t in T])                         ## bus vpltage
@variable(GridModel, θ[i in I, y in Y, d in D, t in T])                         ## bus voltage angle

@variable(GridModel, Ps[s in S,  y in Y, d in D, t in T] >= 0)                  ## active power of solar unit

@variable(GridModel, Pw[w in W,  y in Y, d in D, t in T] >= 0)                  ## active power of wind unit
#@variable(GridModel, Qw[w in W,  y in Y, d in D, t in T])                       ## reactive power of wind unit

@variable(GridModel, Pc[c in C,  y in Y, d in D, t in T])                       ## converter active power
@variable(GridModel, Qc[c in C,  y in Y, d in D, t in T])                       ## converter reactive power

@variable(GridModel, Pl_ac[l in L,  y in Y, d in D, t in T])                    ## active power of a line for AC side
@variable(GridModel, Ql_ac[l in L,  y in Y, d in D, t in T])                    ## reactive power of a line for AC side
@variable(GridModel, Sl_ac[l in L, y in Y, d in D, t in T])                     ## apparent power of line for AC side

@variable(GridModel, Pl_dc[l in L,  y in Y, d in D, t in T])                    ## apparent power of line for DC side

@variable(GridModel, P_e_d[i in I, y in Y, d in D, t in T] >= 0)                ## active load served
@variable(GridModel, Q_e_d[i in I, y in Y, d in D, t in T] >= 0)                ## reactive load served
@constraint(GridModel,P_e_d_max[i in I, y in Y, d in D, t in T], P_e_d[i,y,d,t] <= Pd_max[y]*demandData[string(i)][string(d)][string(t)]["demand"]) ##constraint 1a
@constraint(GridModel,Q_e_d_max[i in I, y in Y, d in D, t in T], Q_e_d[i,y,d,t] <= Pd_max[y]*Pf*demandData[string(i)][string(d)][string(t)]["demand"]) ##constraint 1a

#EV parameters
@variable(GridModel, I_ch[e in E,t in T, d in D, y in Y] ,Bin)                  ## binary variable EV charging
@variable(GridModel, I_chs[ch in CH, y in Y] ,Bin)                              ## binary variable for charging station
@variable(GridModel, I_v2g[e in E,t in T, d in D, y in Y] ,Bin)                 ## binary variable EV V2G
@variable(GridModel, P_ch[e in E,t in T, d in D, y in Y])                       ## EV charging power
@variable(GridModel, P_V2G[e in E,t in T, d in D, y in Y] )                     ## EV V2G power
@variable(GridModel, P_disch[e in E,t in T, d in D, y in Y])                    ## EV discharge power
@variable(GridModel, EV[e in E,t in T, d in D, y in Y])                         ## energy stored in each EV


    ##@variable(GridModel, P_inj_ac[i in I, t in T, d in D, y in Y])                  ## total active power injected for AC side
    ##@variable(GridModel, Q_inj_ac[i in I, t in T, d in D, y in Y])                  ## total rective power injected for AC side

## Objective Function Constraints defined for binary variables

##@cobstraint(GridModel,var_Is[s in S , y in Yn],Is[s,y]>=Is[s, y=!1 ? y-1 : 2)

@constraint(GridModel,var_Iw[w in W , y in Yn],Iw[w,y]>=Iw[w, y-1]) ##constraint 1b
@constraint(GridModel,var_Ic[c in C , y in Yn],Ic[c,y]>=Ic[c, y-1 ]) ##constraint 1c
@constraint(GridModel,var_Il[l in L , y in Yn],Il[l,y]>=Il[l, y-1 ])
@constraint(GridModel,var_Is[s in S , y in Yn],Is[s,y]>=Is[s, y-1 ])
@constraint(GridModel,var_Ich[ch in CH , y in Yn],I_chs[ch,y]>=I_chs[ch, y-1 ])

# Bus voltage and angle constraint
@constraint(GridModel,Var_V[i in I, y in Y, d in D, t in T], Vmin[i] <= V[i,y,d,t] <= Vmax[i])              ##constraint 2a
#@constraint(GridModel,Var_θ[i in I, y in Y, d in D, t in T], -3.14 <= θ[i,y,d,t] <= 3.14)                       ##constraint 2a
    # Converter active power constraint
@constraint(GridModel,Var_Pc_min[c in C, y in Y, d in D, t in T], Pc_min[c]*Ic[c,y] <= Pc[c,y,d,t])         ##constraint 2b
@constraint(GridModel,Var_Pc_max[c in C, y in Y, d in D, t in T], Pc_max[c]*Ic[c,y] >= Pc[c,y,d,t])         ##constraint 2b
    # Converter reactive power constraint
@constraint(GridModel,Var_Qc_min[c in C, y in Y, d in D, t in T], Qc_min[c]*Ic[c,y] <= Qc[c,y,d,t])         ##constraint 2c
@constraint(GridModel,Var_Qc_max[c in C, y in Y, d in D, t in T], Qc_max[c]*Ic[c,y] >= Qc[c,y,d,t])         ##constraint 2c
    # Line apparant power constraint for AC
@constraint(GridModel,Var_Sl_min[l in L, y in Y, d in D, t in T; AC_DC[l] != 1], -SL_MAX[l]*Il[l,y] <= Sl_ac[l,y,d,t])     ##constraint 2e (l=1 to 6)
@constraint(GridModel,Var_Sl_max[l in L, y in Y, d in D, t in T; AC_DC[l] != 1],  SL_MAX[l]*Il[l,y] >= Sl_ac[l,y,d,t])     ##constraint 2e
    # Line apparant power constraint for DC
@constraint(GridModel,Var_Pl_min[l in L, y in Y, d in D, t in T; AC_DC[l] == 1], -SL_MAX[l]*Il[l,y] <= Pl_dc[l,y,d,t])     ##constraint 2f l=6 to 12)
@constraint(GridModel,Var_Pl_max[l in L, y in Y, d in D, t in T; AC_DC[l] == 1],  SL_MAX[l]*Il[l,y] >= Pl_dc[l,y,d,t])     ##constraint 2f
    # Line apparant power constraint for AC
@constraint(GridModel,Var_Sl_ac[l in L, t in T, d in D, y in Y; AC_DC[l] != 1], Sl_ac[l,y,d,t] == Pl_ac[l,y,d,t] + (0.8*Ql_ac[l,y,d,t]))
    ## Renewable Generation Constraints
@constraint(GridModel,Var_Pw[w in W, y in Y, d in D, t in T], Pw[w,y,d,t] <= Pw_cap[w]*(windData[string(d)][string(t)]["wind_curve"])*Iw[w,y]) ##constraint 3a
@constraint(GridModel,Var_Ps[s in S, y in Y, d in D, t in T], Ps[s,y,d,t] <= Ps_cap[s]*(solarData[string(d)][string(t)]["solar_curve"])*Is[s ,y])       ##constraint 3b


## Power Flow Constraints for AC side

@constraint(GridModel, Pl_ac_limit_min[l in L,  y in Y, d in D, t in T; AC_DC[l] != 1], # 1 means the line is a DC one
                    Pl_ac[l,y,d,t] <=
                         -GG[line_from_bus[l],line_to_bus[l]]*(V[line_from_bus[l],y,d,t]-V[line_to_bus[l],y,d,t])
                        + BB[line_from_bus[l],line_to_bus[l]]*(θ[line_from_bus[l],y,d,t]-θ[line_to_bus[l],y,d,t]) + M * (1-Il[l,y])
                        )

@constraint(GridModel, Pl_ac_limit_max[l in L,  y in Y, d in D, t in T; AC_DC[l] != 1],
                    Pl_ac[l,y,d,t] >=
                         -GG[line_from_bus[l],line_to_bus[l]]*(V[line_from_bus[l],y,d,t]-V[line_to_bus[l],y,d,t])
                        + BB[line_from_bus[l],line_to_bus[l]]*(θ[line_from_bus[l],y,d,t]-θ[line_to_bus[l],y,d,t]) - M * (1-Il[l,y])
                        )

@constraint(GridModel, Ql_ac_limit_min[l in L, y in Y, d in D, t in T; AC_DC[l] != 1],
                    Ql_ac[l,y,d,t] <=
                          BB[line_from_bus[l],line_to_bus[l]]*(θ[line_from_bus[l],y,d,t]-θ[line_to_bus[l],y,d,t])
                        + GG[line_from_bus[l],line_to_bus[l]]*(V[line_from_bus[l],y,d,t]-V[line_to_bus[l],y,d,t] + M * (1-Il[l,y]))
                         )

@constraint(GridModel, Ql_ac_limit_max[l in L, y in Y, d in D, t in T; AC_DC[l] != 1],
                     Ql_ac[l,y,d,t] >=
                          BB[line_from_bus[l],line_to_bus[l]]*(θ[line_from_bus[l],y,d,t]-θ[line_to_bus[l],y,d,t])
                        + GG[line_from_bus[l],line_to_bus[l]]*(V[line_from_bus[l],y,d,t]-V[line_to_bus[l],y,d,t] - M * (1-Il[l,y]))
                        )

## Power Flow Constraints for DC side
@constraint(GridModel, p_line_dc1[l in L, y in Y, d in D, t in T; AC_DC[l] == 1], Pl_dc[l,y,d,t] <=
                            (V[line_from_bus[l],y,d,t]-V[line_to_bus[l],y,d,t])/R[l] + M* (1-Il[l,y]) )

@constraint(GridModel, p_line_dc2[l in L, y in Y, d in D, t in T; AC_DC[l] == 1], Pl_dc[l,y,d,t] >=
                            (V[line_from_bus[l],y,d,t]-V[line_to_bus[l],y,d,t])/R[l] - M* (1-Il[l,y]))

## EV Constraints
@constraint(GridModel,Var_Pch_min[e in E, t in T, d in D, y in Y], 0 <= P_ch[e,t,d,y] )
@constraint(GridModel,Var_Pch_max[ch in CH,e in E, t in T, d in D, y in Y], P_ch[e,t,d,y]<= I_ch[e,t,d,y]* Pch_max[e]) #I_e is parametr incidance matrix

@constraint(GridModel,Var_PV2G_min[e in E, t in T, d in D, y in Y],0 <= P_V2G[e,t,d,y] )    ##constraint 6b
@constraint(GridModel,Var_PV2G_max[ch in CH,e in E, t in T, d in D, y in Y],P_V2G[e,t,d,y] <= I_v2g[e,t,d,y]*Pv2g_max[e])

@constraint(GridModel,Var_Ichs[e in E, t in T, d in D, y in Y], I_ch[e,t,d,y] + I_v2g[e,t,d,y]<= sum((Incedene_Dictionary[string(t)][string(e)][string(ch)]["Incidence"])*I_chs[ch,y] for ch in CH))               ##constraint 6c

@constraint(GridModel,Var_Pdisch_min[e in E, t in T, d in D, y in Y], 0 <= P_disch[e,t,d,y] )
@constraint(GridModel,Var_Pdisch_max[e in E, t in T, d in D, y in Y], P_disch[e,t,d,y] == Idisch[e,t] * Ptravel[d][t] * 0.3/1000) # from kWh to MWh --> EV discharge should be in p.u

@constraint(GridModel,Energy_Stored_EV_min[e in E, t in T, d in D, y in Y], EV_min[e] <= EV[e,t,d,y] )
@constraint(GridModel,Energy_Stored_EV_max[e in E, t in T, d in D, y in Y], EV_max[e] >= EV[e,t,d,y])


@constraint(GridModel,EnergyStoredPerEV[e in E, t in Tn, d in D, y in Y], EV[e,t,d,y] ==
                EV[e, t!=1 ? t-1 : 24 ,d,y] + P_ch[e,t,d,y] - P_disch[e,t,d,y] - P_V2G[e,t,d,y])           ##constraint 6d


## Nodal Balance Constraints

@constraint(GridModel, Balance_ac_active[i in I, y in Y, d in D, t in T; AC_DC[i] != 1],
                            P_e_d[i,y,d,t] ==
                            sum((P_V2G[e,t,d,y] - P_ch[e,t,d,y]) for e in E if EV_bus[e]==i)
                            +sum(Ps[s,y,d,t] for s in S if solar_bus[s]==i)
                            +sum(Pw[w,y,d,t] for w in W if wind_bus[w]==i)
                            +sum(Pc[c,y,d,t] for c in C if c_bus_ac[c]==i)
                            +sum(Pl_ac[l,y,d,t] for l in L if line_to_bus[l] == i)
                            -sum(Pl_ac[l,y,d,t] for l in L if line_from_bus[l] == i))       ##(v2g -ch)         ##constraint 5a

@constraint(GridModel, Balance_ac_reactive[i in I,  y in Y, d in D, t in T; AC_DC[i] != 1],
                            Q_e_d[i,y,d,t] ==
                            +sum(Qc[c,y,d,t] for c in C if c_bus_ac[c]==i)
                            +sum(Ql_ac[l,y,d,t] for l in L if line_to_bus[l] == i)
                            -sum(Ql_ac[l,y,d,t] for l in L if line_from_bus[l] == i))

@constraint(GridModel, Balance_dc_active[e in E, i in I,  y in Y, d in D, t in T; AC_DC[i] == 1],
                            P_e_d[i,y,d,t] ==
                            sum((P_V2G[e,t,d,y] - P_ch[e,t,d,y]) for e in E if EV_bus[e]==i)
                            +sum(Ps[s,y,d,t] for s in S if solar_bus[s]==i)
                            +sum(Pw[w,y,d,t] for w in W if wind_bus[w]==i)
                            -sum(Pc[c,y,d,t] for c in C if c_bus_dc[c]==i)
                            +sum(Pl_dc[l,y,d,t] for l in L if line_to_bus[l] == i)
                            -sum(Pl_dc[l,y,d,t] for l in L if line_from_bus[l] == i) )                                                                                                              ##constraint 5b

## Objective function
@objective(GridModel, Min, (sum(cost_solar/((1+r)^(y-1)) * (Is[s ,y] - Is[s, y-1]) for s in S, y in Yn) + sum(cost_solar * (Is[s ,1]) for s in S )
                            +sum(cost_wind/((1+r)^(y-1)) * (Iw[w ,y] - Iw[w, y-1]) for w in W, y in Yn) + sum(cost_wind * (Iw[w ,1]) for w in W )
                            +sum(cost_converter/((1+r)^(y-1)) * (Ic[c ,y] - Ic[c, y-1]) for c in C, y in Yn) + sum(cost_converter * (Ic[c ,1]) for c in C )
                            +sum(cost_line/((1+r)^(y-1)) * (Il[l ,y] - Il[l, y-1]) for l in L, y in Yn) + sum(cost_line * (Il[l ,1]) for l in L )
                            +sum(cost_ChargingStation/((1+r)^(y-1)) * (I_chs[ch ,y] - I_chs[ch, y-1]) for ch in CH, y in Yn) + sum(cost_ChargingStation * (I_chs[ch ,1]) for ch in CH )
                            +K*sum(Pd_max[y]*demandData[string(i)][string(d)][string(t)]["demand"] - P_e_d[i,y,d,t] for i in I,  y in Y, d in D, t in T)))

optimize!(GridModel)
obj = objective_value(GridModel)

Iw_sol = JuMP.value.(Iw)
Is_sol = JuMP.value.(Is)
Il_sol = JuMP.value.(Il)
Ic_sol = JuMP.value.(Ic)
I_chs_sol = JuMP.value.(I_chs)


#=
sum((JuMP.value(Ps[s,y,d,t]) - (Ps_cap[s]*(solarData[string(d)][string(t)]["solar_curve"])*JuMP.value(Is[s,y]))) for s in S,  y in 1 , d in 1, t in T)
sum(JuMP.value(Pw[w,y,d,t]) - (Pw_cap[w]*(windData[string(d)][string(t)]["wind_curve"])*JuMP.value(Iw[w,y])) for w in W,  y in 1 , d in 1, t in T)
sum((Ps_cap[s]*(solarData[string(d)][string(t)]["solar_curve"])*JuMP.value(Is[s,y])) for s in S,  y in 1 , d in 1, t in T)
sum(JuMP.value(Ps[s,1,1,12]) for s in S)+ sum(JuMP.value(Pw[w,1,1,12]) for w in W)+ sum(JuMP.value(P_V2G[e,12,1,1]) - JuMP.value(P_ch[e,12,1,1]) for e in E)
sum(JuMP.value(P_e_d[i,1,1,12]) for i in I)
sum(Pd_max[y]*demandData[string(i)][string(d)][string(t)]["demand"] - JuMP.value(P_e_d[i,y,d,t]) for i in I,  y in 10, d in 1:7, t in 1:24)
penalty = K*sum(Pd_max[y]*demandData[string(i)][string(d)][string(t)]["demand"] - JuMP.value(P_e_d[i,y,d,t]) for i in I,  y in Y, d in D, t in T)
display(sum(cost_solar/((1+r)^(y-1)) * (Is[s ,y] - Is[s, y-1]) for s in S, y in Yn) + sum(cost_solar * (Is[s ,1]) for s in S ))
Ped=sum(Pd_max[y]*demandData[string(i)][string(d)][string(t)]["demand"]  for i in I,  y in Y, d in D, t in T)
=#
