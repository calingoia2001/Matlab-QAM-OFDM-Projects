function[rounded_values]=Roundtospecific(values,specific)
rounded_values=values;
    for d = 1: length(values)
        v=abs(specific-values(d));
        [~,index]=min(v);
        rounded_values(d)=specific(index);
        
    end
end