function sum_var=nnl_serial1(term,start_N,sum_N,step)
    sum_var=0;
    for i=start_N:step:sum_N
        sum_var=sum_var+single(term(i));
    end
end
