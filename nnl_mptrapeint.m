function int_res=nnl_mptrapeint(span,points_n,f_expr)
    points=linspace(span(1),span(2),points_n);
    int_res=(span(2)-span(1))/(2*points_n)*(f_expr(span(1))+f_expr(span(2))+2*sum(arrayfun(f_expr,points(1:points_n-1))));
end
