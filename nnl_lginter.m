function interpfunc=nnl_lginter(span,points_n,f_expr,var_vec)
    if length(span)~=2
        error('SPAN PARAMETER MEUST BE LIKE [A,B]')
        if span(1)>span(2)
            error('SPAN [A,B],A MUST BE LESS THAN B')
        end
    end
    
    points=linspace(span(1),span(2),points_n);

    if isa(f_expr,'function_handle')
        f_points=arrayfun(f_expr,points);
    else
        f_points=f_expr;
    end

    f_items={points_n};
    j=1;
    for i=points
        x1vec=points(1:end ~=j);
        x2vec=i-x1vec;
        f_items{j,1}=@(x)prod(x-x1vec)/(prod(x2vec));
        j=j+1;
    end

    interpfunc=zeros(1,length(var_vec));
    for i=1:points_n
        interpfunc=interpfunc+arrayfun(f_items{i},var_vec)*f_points(i);
    end
end
