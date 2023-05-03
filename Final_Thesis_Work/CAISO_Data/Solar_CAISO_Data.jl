import.pkg
Pkg.add("HTTP")
Pkg.add("ZipFile")
Pkg.add("InfoZIP")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("LightXML")
using HTTP
using ZipFile
using InfoZIP
using CSV, DataFrames, LightXML
using Dates
loc=pwd()

# Get the path of the current directory (where my_script.jl is located)
current_dir = @__DIR__

# The type of data that we want is the Market/Process and it can be one of these options: DAM (DAM prices), ACTUAL,RTM,2DA,7DA,RTDP,RTD (Real Time Dispatch),HASP

##solar MW data for NP15, SP15, ZP26
Solar_Years_NP15 = DataFrame(INTERVALSTARTTIME_GMT= String31[], INTERVALENDTIME_GMT=String31[],MW =Float64[] )
Solar_Years_SP15 = DataFrame(INTERVALSTARTTIME_GMT= String31[], INTERVALENDTIME_GMT=String31[],MW =Float64[] )
Solar_Years_ZP26 = DataFrame(INTERVALSTARTTIME_GMT= String31[], INTERVALENDTIME_GMT=String31[],MW =Float64[] )

## Start date jan 2019
start=Dates.DateTime(2020, 12, 31, 00, 00)
HoursOfDays =24
BunchOfDays = 30
DaysOfYear = 365
NumberOfYears = 2

TotalNumberOfDays = NumberOfYears*DaysOfYear
NumberOfIterrations = floor((TotalNumberOfDays/BunchOfDays)+1)

Days = 1

## download 37*30 days data and store it in 3 different arrays
for iter in 1:NumberOfIterrations
    "solar"
    beg = start + Day(Days)
    beg_time = Dates.format(beg,"yyyymmddTHH:MM")
    en = beg + Day(BunchOfDays)
    end_time = Dates.format(en,"yyyymmddTHH:MM")

    QN = "SLD_REN_FCST"
    MKT = "ACTUAL"

    url_renewable_act="http://oasis.caiso.com/oasisapi/SingleZip?resultformat=6&queryname=$(QN)&version=1&market_run_id=$(MKT)&startdatetime=$(string(beg_time))-0000&enddatetime=$(string(end_time))-0000"
    name_down=HTTP.download(url_renewable_act, loc)

    sleep(5)

    InfoZIP.unzip(name_down, loc)

    println(chop(name_down,head=0,tail=4))

    csv_reader = CSV.File("$(chop(name_down,head=0,tail=4))CSV")

    df_renewable = DataFrame(csv_reader)

    #println(df_renewable)

##  NP15 data
    solar_bunchOfDays_NP15_AllColumns = filter(row -> (row["TRADING_HUB"] == "NP15") & (row["RENEWABLE_TYPE"] == "Solar"), df_renewable)
    solar_bunchOfDays_NP15 = solar_bunchOfDays_NP15_AllColumns[:, [4,5,12]]
    solar_bunchOfDays_NP15_Sorted = sort(solar_bunchOfDays_NP15, order(:INTERVALENDTIME_GMT, ))

## SP15 data
    solar_bunchOfDays_SP15_AllColumns = filter(row -> (row["TRADING_HUB"] == "SP15") & (row["RENEWABLE_TYPE"] == "Solar"), df_renewable)
    solar_bunchOfDays_SP15 = solar_bunchOfDays_SP15_AllColumns[:, [4,5,12]]
    solar_bunchOfDays_SP15_Sorted = sort(solar_bunchOfDays_SP15, order(:INTERVALENDTIME_GMT, ))

##  ZP26 data
    solar_bunchOfDays_ZP26_AllColumns = filter(row -> (row["TRADING_HUB"] == "ZP26") & (row["RENEWABLE_TYPE"] == "Solar"), df_renewable)
    solar_bunchOfDays_ZP26 = solar_bunchOfDays_ZP26_AllColumns[:, [4,5,12]]
    solar_bunchOfDays_ZP26_Sorted = sort(solar_bunchOfDays_ZP26, order(:INTERVALENDTIME_GMT, ))

## finding and replacing missing hour data
    for i in 1:((BunchOfDays * HoursOfDays) - 1)

        pn = DateTime("$(chop(solar_bunchOfDays_NP15_Sorted[i,1] ,head=0,tail=6))")
        nn = DateTime("$(chop(solar_bunchOfDays_NP15_Sorted[i+1,1] ,head=0,tail=6))")
        if ((Hour)(nn-pn) != Hour(1))
            mean = (solar_bunchOfDays_NP15_Sorted[i,3] + solar_bunchOfDays_NP15_Sorted[i+1,3])/2
            insert!.(eachcol(solar_bunchOfDays_NP15_Sorted), i+1, [replace("$(Dates.format(pn+Hour(1), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"), replace("$(Dates.format(pn+Hour(2), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"),mean])
        end

        ps = DateTime("$(chop(solar_bunchOfDays_SP15_Sorted[i,1] ,head=0,tail=6))")
        ns = DateTime("$(chop(solar_bunchOfDays_SP15_Sorted[i+1,1] ,head=0,tail=6))")
        if ((Hour)(ns-ps) != Hour(1))
            mean = (solar_bunchOfDays_SP15_Sorted[i,3] + solar_bunchOfDays_SP15_Sorted[i+1,3])/2
            insert!.(eachcol(solar_bunchOfDays_SP15_Sorted), i+1, [replace("$(Dates.format(ps+Hour(1), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"), replace("$(Dates.format(ps+Hour(2), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"),mean])
        end

        pz = DateTime("$(chop(solar_bunchOfDays_ZP26_Sorted[i,1] ,head=0,tail=6))")
        nz = DateTime("$(chop(solar_bunchOfDays_ZP26_Sorted[i+1,1] ,head=0,tail=6))")
        if ((Hour)(nz-pz) != Hour(1))
            mean = (solar_bunchOfDays_ZP26_Sorted[i,3] + solar_bunchOfDays_ZP26_Sorted[i+1,3])/2
            insert!.(eachcol(solar_bunchOfDays_ZP26_Sorted), i+1, [replace("$(Dates.format(pz+Hour(1), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"), replace("$(Dates.format(pz+Hour(2), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"),mean])
        end

    end

## appending 30 days daata on each iteration
    append!(Solar_Years_NP15, solar_bunchOfDays_NP15_Sorted)
    append!(Solar_Years_SP15, solar_bunchOfDays_SP15_Sorted)
    append!(Solar_Years_ZP26, solar_bunchOfDays_ZP26_Sorted)

    Days = Days + BunchOfDays;

end

## sort the gathered data
Final_Solar_Years_NP15 = sort(Solar_Years_NP15, order(:INTERVALENDTIME_GMT, ))
println(Final_Solar_Years_NP15)
println(Final_Solar_Years_NP15[1417,:])

csv_file_path_solar_NP15= joinpath(current_dir, "Demand_Datafiles", "Final_Solar_Years_NP15.csv")
CSV.write(csv_file_path_solar_NP15, Final_Solar_Years_NP15)

## remove NP15 data for 29th feb 2020
np = 1417
for i in 1 : 24
    deleteat!(Final_Solar_Years_NP15, np)
end
println(Final_Solar_Years_NP15[1417,:])
final=first(Final_Solar_Years_NP15,1427)
println(final)

Final_Solar_Years_SP15 = sort(Solar_Years_SP15, order(:INTERVALENDTIME_GMT, ))
println(Final_Solar_Years_SP15)
sp = 1417
for i in 1 : 24
    deleteat!(Final_Solar_Years_SP15, sp)
end
println(Final_Solar_Years_SP15[1417,:])

csv_file_path_solar_SP15= joinpath(current_dir, "Demand_Datafiles", "Final_Solar_Years_SP15.csv")
CSV.write(csv_file_path_solar_SP15, Final_Solar_Years_SP15)

Final_Solar_Years_ZP26 = sort(Solar_Years_ZP26, order(:INTERVALENDTIME_GMT, ))
println(Final_Solar_Years_ZP26)
zp = 1417
for i in 1 : 24
    deleteat!(Final_Solar_Years_ZP26, zp)
end
println(Final_Solar_Years_ZP26[1417,:])

csv_file_path_solar_ZP26 = joinpath(current_dir, "Demand_Datafiles", "Final_Solar_Years_ZP26.csv")
CSV.write(csv_file_path_solar_ZP26, Final_Solar_Years_ZP26)

## removing negative values of NP15
n=nrow(Final_Solar_Years_NP15)
for i in 1 : n
    if (Final_Solar_Years_NP15[i,3] < 0)
        Final_Solar_Years_NP15[i,3] = 0.00
    end
end
println(Final_Solar_Years_NP15)

j = 1
## first year data matrix for NP15
days = 1
Solar_NP15_1 = Array{Float64, 2}(undef, 365, 24)
for i in 1:365

        Solar_NP15_1[days,:]= Final_Solar_Years_NP15[j:j+23,3]
        j=j+24
        days = days+1

end
display(Solar_NP15_1)
Solar_df_NP15_1=DataFrame(Solar_NP15_1, :auto)

csv_file_path_Solar_NP15_1 = joinpath(current_dir, "Demand_Datafiles", "Solar_df_NP15_1.csv")
CSV.write(csv_file_path_Solar_NP15_1, Solar_df_NP15_1)
## removing negative values of SP15
n=nrow(Final_Solar_Years_SP15)
for i in 1 : n
    if (Final_Solar_Years_SP15[i,3] < 0)
        Final_Solar_Years_SP15[i,3] = 0.00
    end
end
println(Final_Solar_Years_SP15)
j = 1

## first year data matrix for SP15
days = 1
Solar_SP15_1 = Array{Float64, 2}(undef, 365, 24)
for i in 1:365

        Solar_SP15_1[days,:]= Final_Solar_Years_SP15[j:j+23,3]
        j=j+24
        days = days+1
end
display(Solar_SP15_1)
Solar_df_SP15_1=DataFrame(Solar_SP15_1, :auto)

csv_file_path_Solar_SP15_1 = joinpath(current_dir, "Demand_Datafiles", "Solar_df_SP15_1.csv")
CSV.write(csv_file_path_Solar_SP15_1, Solar_df_SP15_1)
## removing negative values of ZP26
n=nrow(Final_Solar_Years_ZP26)
for i in 1 : n
    if (Final_Solar_Years_ZP26[i,3] < 0)
        Final_Solar_Years_ZP26[i,3] = 0.00
    end
end
println(Final_Solar_Years_ZP26)
j = 1

## First year data matrix for ZP26
days = 1
Solar_ZP26_1 = Array{Float64, 2}(undef, 365, 24)
for i in 1:365.

        Solar_ZP26_1[days,:]= Final_Solar_Years_ZP26[j:j+23,3]
        j=j+24
        days = days+1
end
display(Solar_ZP26_1)
Solar_df_ZP26_1=DataFrame(Solar_ZP26_1, :auto)
csv_file_path_Solar_ZP26_1 = joinpath(current_dir, "Demand_Datafiles", "Solar_df_ZP26_1.csv")
CSV.write(csv_file_path_Solar_ZP26_1, Solar_df_ZP26_1)
## First Year SUMMATION OF NP,SP,ZP
SOLAR_NP_SP_ZP_1=Solar_NP15_1 + Solar_SP15_1 + Solar_ZP26_1
SOLAR_NP_SP_ZP_1_df = DataFrame(SOLAR_NP_SP_ZP_1, :auto)
println(SOLAR_NP_SP_ZP_1_df)

csv_file_path_SOLAR_NP_SP_ZP_1_df = joinpath(current_dir, "Demand_Datafiles", "SOLAR_NP_SP_ZP_1_df.csv")
CSV.write(csv_file_path_SOLAR_NP_SP_ZP_1_df, SOLAR_NP_SP_ZP_1_df)
