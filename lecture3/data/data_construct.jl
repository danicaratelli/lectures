function data_construct(Mnems,Tran,start_date,vint,Nseries);
    Data_ = repmat([NaN],1000,Nseries);
    dts = [];
    dates = [];
    quarts =  Integer[];
    q = 1;
    for k=1:Nseries;
        x = get_data(f,Mnems[k],vintage_dates=vint,units=Tran[k]);
        dts = x.df[:date];
        data = x.df[:value];
        st = find(dts.==Date(start_date));
        if isempty(st);
            dates_add = (Date(Dates.year(Date(start_date))):
                         Date(dts[1]));
            dates_add= (recur(dates_add) do x
                        Dates.dayofmonth(x) == 1
                        end);
            data = [repmat([NaN],length(dates_add)-1,1);data];
            dts = [dates_add; dts[2:end]];
        end
        st = find(dts.==Date(start_date));
        if x.freq=="Monthly";
            dates = dts[st[1]:end];
        elseif x.freq=="Quarterly";
            push!(quarts,k);
            q = q+1;
        end
        data = data[st[1]:end];
        if x.freq=="Quarterly";
            data = map(x->[NaN;NaN;x],data);
            data = vcat([data[x] for x in 1:length(data)]...);
        end
        Data_[1:size(data,1),k] = reshape(data,length(data),1);
    end
    #moving quarterly series at the end
    Data_ = [Data_[:,setdiff(1:Nseries,quarts)] Data_[:,quarts]];
    return Data(dates,Data_,quarts);
end
