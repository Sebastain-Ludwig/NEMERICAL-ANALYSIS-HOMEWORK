function int_res=nnl_mpsimpsonint(span,points_n,f_expr)
    points=linspace(span(1),span(2),points_n*2);
    int_res=(span(2)-span(1))/(6*points_n)*(f_expr(span(1))+f_expr(span(2))+4*sum(arrayfun(f_expr,points(1:2:end-1)))+2*sum(arrayfun(f_expr,points(2:2:end))));
end
