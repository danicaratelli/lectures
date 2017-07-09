## Packages used
using Dierckx, FredData, Dates, DataFrames

##Optionsf
FRED_API_key = "4f9f790c74d8b1e1f808400bc343fbe3"; #St. Louis Fed Alfred jey

start_date = "2000-01-01";

date0 = "2017-02-28"; date0 = Date(date0);
date1 = "2017-03-31"; date1 = Date(date1);
date2 = "2017-04-30"; date2 = Date(date2);
date3 = "2017-05-31"; date3 = Date(date3);

type Data
    dates::Array{Date,1};
    Data:: Array{Float64,2};
    quarts:: Array{Int64,1};
end


##Setting prerequisites

#useful functions
include("data_construct.jl");
include("menmonics.jl");
include("transformation.jl");
include("series.jl");

#Fred pulling
f = Fred();

Nseries = length(Mnems);
#data 1
data0 = data_construct(Mnems,Transformation,start_date,date0,Nseries);
data1 = data_construct(Mnems,Transformation,start_date,date1,Nseries);
data2 = data_construct(Mnems,Transformation,start_date,date2,Nseries);
data3 = data_construct(Mnems,Transformation,start_date,date3,Nseries);
