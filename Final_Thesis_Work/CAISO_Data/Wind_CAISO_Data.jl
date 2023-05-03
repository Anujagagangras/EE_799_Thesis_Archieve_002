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
# The type of data that we want is the Market/Process and it can be one of these options: DAM (DAM prices), ACTUAL,RTM,2DA,7DA,RTDP,RTD (Real Time Dispatch),HASP

# Get the path of the current directory (where my_script.jl is located)
current_dir = @__DIR__
##solar MW data for NP15, SP15, ZP26
Wind_NP15 = DataFrame(INTERVALSTARTTIME_GMT= String31[], INTERVALENDTIME_GMT=String31[],MW =Float64[] )
Wind_SP15 = DataFrame(INTERVALSTARTTIME_GMT= String31[], INTERVALENDTIME_GMT=String31[],MW =Float64[] )

## Start date jan 2019
Start=Dates.DateTime(2020, 12, 31, 00, 00)
HoursOfDays =24
BunchOfDays = 30
DaysOfYear = 365
NumberOfYears = 2
TotalNumberOfDays = NumberOfYears*DaysOfYear
NumberOfIterrations = floor((TotalNumberOfDays/BunchOfDays)+1)

Days = 1
## download 37*30 days data and store it in 3 different arrays
for iter in 1:NumberOfIterrations
    "Wind"
    beg = Start + Day(Days)
    beg_time = Dates.format(beg,"yyyymmddTHH:MM")
    en = beg + Day(BunchOfDays)
    end_time = Dates.format(en,"yyyymmddTHH:MM")

    QN = "SLD_REN_FCST"
    MKT = "ACTUAL"

    url_renewable_act="http://oasis.caiso.com/oasisapi/SingleZip?resultformat=6&queryname=$(QN)&version=1&market_run_id=$(MKT)&startdatetime=$(string(beg_time))-0000&enddatetime=$(string(end_time))-0000"
    name_down=HTTP.download(url_renewable_act, loc)

    sleep(3)

    InfoZIP.unzip(name_down, loc)

    println(chop(name_down,head=0,tail=4))

    csv_reader = CSV.File("$(chop(name_down,head=0,tail=4))CSV")

    wind_data = DataFrame(csv_reader)

##  NP15 data
    wind_bunchOfDays_NP15_AllColumns = filter(row -> (row["TRADING_HUB"] == "NP15") & (row["RENEWABLE_TYPE"] == "Wind"), wind_data)
    wind_bunchOfDays_NP15 = wind_bunchOfDays_NP15_AllColumns[:, [4,5,12]]
    wind_bunchOfDays_NP15_Sorted = sort(wind_bunchOfDays_NP15, order(:INTERVALENDTIME_GMT, ))

## SP15 data
    wind_bunchOfDays_SP15_AllColumns = filter(row -> (row["TRADING_HUB"] == "SP15") & (row["RENEWABLE_TYPE"] == "Wind"), wind_data)
    wind_bunchOfDays_SP15 = wind_bunchOfDays_SP15_AllColumns[:, [4,5,12]]
    wind_bunchOfDays_SP15_Sorted = sort(wind_bunchOfDays_SP15, order(:INTERVALENDTIME_GMT, ))

## finding and replacing missing hour data
    for i in 1:((BunchOfDays * HoursOfDays) - 1)

        pn = DateTime("$(chop(wind_bunchOfDays_NP15_Sorted[i,1] ,head=0,tail=6))")
        nn = DateTime("$(chop(wind_bunchOfDays_NP15_Sorted[i+1,1] ,head=0,tail=6))")
        if ((Hour)(nn-pn) != Hour(1))
            mean = (wind_bunchOfDays_NP15_Sorted[i,3] + wind_bunchOfDays_NP15_Sorted[i+1,3])/2
            insert!.(eachcol(wind_bunchOfDays_NP15_Sorted), i+1, [replace("$(Dates.format(pn+Hour(1), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"), replace("$(Dates.format(pn+Hour(2), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"),mean])
        end

        ps = DateTime("$(chop(wind_bunchOfDays_SP15_Sorted[i,1] ,head=0,tail=6))")
        ns = DateTime("$(chop(wind_bunchOfDays_SP15_Sorted[i+1,1] ,head=0,tail=6))")
        if ((Hour)(ns-ps) != Hour(1))
            mean = (wind_bunchOfDays_SP15_Sorted[i,3] + wind_bunchOfDays_SP15_Sorted[i+1,3])/2
            insert!.(eachcol(wind_bunchOfDays_SP15_Sorted), i+1, [replace("$(Dates.format(ps+Hour(1), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"), replace("$(Dates.format(ps+Hour(2), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"),mean])
        end
    end

## appending 30 days daata on each iteration
    append!(Wind_NP15, wind_bunchOfDays_NP15_Sorted)
    append!(Wind_SP15, wind_bunchOfDays_SP15_Sorted)

    Days = Days + BunchOfDays;

end
## sort the gathered data
Final_Wind_NP15 = sort(Wind_NP15, order(:INTERVALENDTIME_GMT, ))
println(Final_Wind_NP15)
println(Final_Wind_NP15[1417,:])

## remove NP15 data for 29th feb 2020
np = 1417
for i in 1 : 24
    deleteat!(Final_Wind_NP15, np)
end
println(Final_Wind_NP15[1417,:])
final=first(Final_Wind_NP15,1427)
println(final)

Final_Wind_SP15 = sort(Wind_SP15, order(:INTERVALENDTIME_GMT, ))
println(Final_Wind_SP15)
sp = 1417
for i in 1 : 24
    deleteat!(Final_Wind_SP15, sp)
end
println(Final_Wind_SP15[1417,:])
## write sorted data to csv file
csv_file_path_Wind_NP15 = joinpath(current_dir, "Demand_Datafiles", "Final_Wind_NP15.csv")
# Write the CSV file to the specified path
CSV.write(csv_file_path_Wind_NP15, Final_Wind_NP15)

csv_file_path_Wind_SP15 = joinpath(current_dir, "Demand_Datafiles", "Final_Wind_SP15.csv")
CSV.write(csv_file_path_Wind_SP15, Final_Wind_SP15)

## removing negative values of NP15
n=nrow(Final_Wind_NP15)
for i in 1 : n
    if (Final_Wind_NP15[i,3] < 0)
        Final_Wind_NP15[i,3] = 0.00
    end
end
println(Final_Wind_NP15)

j = 1
## first year data matrix for NP15
days = 1
FNP15_1 = Array{Float64, 2}(undef, 365, 24)
for i in 1:365

        FNP15_1[days,:]= Final_Wind_NP15[j:j+23,3]
        j=j+24
        days = days+1
end
display(FNP15_1)
Wind_NP15_1=DataFrame(FNP15_1, :auto)

## removing negative values of SP15
n=nrow(Final_Wind_SP15)
for i in 1 : n
    if (Final_Wind_SP15[i,3] < 0)
        Final_Wind_SP15[i,3] = 0.00
    end
end
println(Final_Wind_SP15)
j = 1
## first year data matrix for SP15
days = 1
FSP15_1 = Array{Float64, 2}(undef, 365, 24)
for i in 1:365

        FSP15_1[days,:]= Final_Wind_SP15[j:j+23,3]
        j=j+24
        days = days+1
end
display(FSP15_1)
Wind_SP15_1=DataFrame(FSP15_1, :auto)

## First Year SUMMATION OF NP,SP wind data
WIND_NP_SP=FNP15_1 + FSP15_1
Wind_NP_SP = DataFrame(WIND_NP_SP, :auto)
println(Wind_NP_SP)

csv_file_path_Wind_NP_SP_1 = joinpath(current_dir, "Demand_Datafiles", "Wind_NP_SP_1.csv")
CSV.write(csv_file_path_Wind_NP_SP_1, Wind_NP_SP_1)
