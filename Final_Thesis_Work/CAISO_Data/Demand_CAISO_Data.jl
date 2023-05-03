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

## Demand DATA for 2020
Demand_Years = DataFrame(INTERVALSTARTTIME_GMT= String31[], INTERVALENDTIME_GMT=String31[],MW =Float64[] )

## Start date jan 2019
Start=Dates.DateTime(2019, 12, 31, 00, 00)
HoursOfDays =24
BunchOfDays = 30
DaysOfYear = 365
NumberOfYears = 2

TotalNumberOfDays = NumberOfYears*DaysOfYear
NumberOfIterrations = floor((TotalNumberOfDays/BunchOfDays)+1)
Days = 1
## download 37*30 days data and store it in array
for iter in 1:NumberOfIterrations
    "Demand"
    beg = Start + Day(Days)
    beg_time=Dates.format(beg,"yyyymmddTHH:MM")
    en=beg + Day(BunchOfDays)
    end_time=Dates.format(en,"yyyymmddTHH:MM")
    QN_FCST,MKT_RUN="SLD_FCST","ACTUAL" #The forecast MW demand for each hour of the Operating Day, posted in the morning the day before the Operating Day, before the markets run;
    url_FCST="http://oasis.caiso.com/oasisapi/SingleZip?resultformat=6&queryname=$(QN_FCST)&market_run_id=$(MKT_RUN)&startdatetime=$(string(beg_time))-0000&enddatetime=$(string(end_time))-0000&version=1"
    name_down=HTTP.download(url_FCST, loc)
    # http://oasis.caiso.com/oasisapi/SingleZip?resultformat=6&queryname=SLD_FCST&version=1&market_run_id=ACTUAL&startdatetime=20220221T08:00-0000&enddatetime=20220322T07:00-0000
    InfoZIP.unzip(name_down, loc)
    sleep(6)


    InfoZIP.unzip(name_down, loc)

    println(chop(name_down,head=0,tail=4))

    csv_reader = CSV.File("$(chop(name_down,head=0,tail=4))CSV")

    df_demand = DataFrame(csv_reader)

    demand_bunchOfDays_AllColumns=filter(row -> (row["TAC_AREA_NAME"] == "CA ISO-TAC"), df_demand)
    demand_bunchOfDays=demand_bunchOfDays_AllColumns[:, [1,2,12]]
    demand_bunchOfDays_sorted = sort(demand_bunchOfDays, order(:INTERVALENDTIME_GMT, ))
    println(demand_bunchOfDays_sorted)
    r=nrow(demand_bunchOfDays_sorted)

## finding and replacing missing hour data
    for i in 1:((BunchOfDays * HoursOfDays) - 1)

        pn = DateTime("$(chop(demand_bunchOfDays_sorted[i,1] ,head=0,tail=6))")
        nn = DateTime("$(chop(demand_bunchOfDays_sorted[i+1,1] ,head=0,tail=6))")
        if ((Hour)(nn-pn) != Hour(1))
            mean = (demand_bunchOfDays_sorted[i,3] + demand_bunchOfDays_sorted[i+1,3])/2
            insert!.(eachcol(demand_bunchOfDays_sorted), i+1, [replace("$(Dates.format(pn+Hour(1), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"), replace("$(Dates.format(pn+Hour(2), "yyyy-mm-dd HH:MM:SS"))-00:00", " "=>"T"),mean])
        end
    end

## appending 30 days data on each iteration
    append!(Demand_Years, demand_bunchOfDays_sorted)
    Days = Days + BunchOfDays;
end

##
println(Demand_Years)
Demand_Years_sorted = sort(Demand_Years, order(:INTERVALENDTIME_GMT, ))
println(Demand_Years_sorted)
csv_file_path_sorted = joinpath(current_dir, "Demand_Datafiles", "Demand_Years_sorted.csv")
CSV.write(csv_file_path_sorted, Demand_Years_sorted)
## remove NP15 data for 29th feb 2020
println(Demand_Years_sorted[1417,:])
np = 1417
for i in 1 : 24
    deleteat!(Demand_Years_sorted, np)
end
println(Demand_Years_sorted[1417,:])
final=first(Demand_Years_sorted,1418)
println(final)
j = 1

##Demand Data for first year
days = 1
demand_1 = Array{Float64, 2}(undef, 365, 24)
for i in 1:365.

        demand_1[days,:]= Demand_Years_sorted[j:j+23,3]
        j=j+24
        days = days+1

end
##Demand Data for second year
days = 1
demand_2 = Array{Float64, 2}(undef, 365, 24)
for i in 1:365

        demand_2[days,:]= Demand_Years_sorted[j:j+23,3]
        j=j+24
        days = days+1

end
display(demand_1)
display(demand_2)
Final_Demand_1=DataFrame(demand_1, :auto)
Final_Demand_2=DataFrame(demand_2, :auto)

# Construct the path to the CSV file relative to the current directory
csv_file_path_d1 = joinpath(current_dir, "Demand_Datafiles", "Final_Demand_1.csv")
# Construct the path to the CSV file relative to the current directory
csv_file_path_d2 = joinpath(current_dir, "Demand_Datafiles", "Final_Demand_2.csv")
# Write the CSV file to the specified path
CSV.write(csv_file_path_d1, Final_Demand_1)
CSV.write(csv_file_path_d2, Final_Demand_2)
