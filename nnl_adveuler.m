function var_traj=nnl_adveuler(var_init,steph,span,f_expr)
    if ~isa(f_expr,'function_handle')
        error('THE INPUT TYPE IS NOT FUNCTION_HANDLE<REQUIRED>');
    end
    
    if length(span)~=2
        error('SPAN LENGTH MUST BE 2')
        if span(1)>span(2)
            error('SPAN:[A,B],A MUST BE LESS THAN B');
        end
    end

    k1=@(t,u,h,his)f_expr(t,u);
    k2=@(t,u,h,his)f_expr(t+h,u+h*his);

    var_traj=nnl_desolve(var_init,steph,span,2,k1,k2);
end
