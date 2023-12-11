function cal_var=hornercal(coevec,var)
%%%HORNER MATHED

%%RECUSIVE FUNCTION
    if isempty(coevec)
        cal_var=1;
    else
        cal_var=coevec(1)+var*hornercal(coevec(2:end),var);
    end
end