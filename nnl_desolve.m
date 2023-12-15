function var_traj=nnl_desolve(var_init,steph,span,f_n,varargin)
    if nargin-4~=f_n
        error('THE DELARE U AND F NUMBER IS NOT EQUAL WITH INPUT');
    else
        for expr=varargin
            if ~isa(expr{1},'function_handle')
                error('THE INPUT TYPE IS NOT FUNCTION_HANDLE<REQUIRED>');
            end
        end
    end

    var_his=repmat(var_init,[1,nargin]);
    var_traj=var_his;
    stepn=0;
    while(span(1)+stepn*steph<span(2))
        hist_f=0;
        sum_f=0;
        for i=1:f_n
            hist_f=varargin{i}(span(1)+steph*stepn,var_traj(end),steph,hist_f);
            sum_f=sum_f+hist_f;
        end
        var_traj=[var_traj,steph*sum_f+var_traj(end)];
        stepn=stepn+1;
    end
end
