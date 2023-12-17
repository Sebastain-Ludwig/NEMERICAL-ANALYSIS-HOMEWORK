function [root,stepn]=nnl_newton_nle(nle_function,nle_derive,var_0,tor_err)
    if abs(nle_derive(var_0))<tor_err/2
        error('THE DERIVE FUNCTION CONT USE NEWTON MATHED,ERROR ON POINT X= %s',var_0);
    end

    var_next=var_0-nle_function(var_0)/nle_derive(var_0);
    
    if abs(nle_function(var_0))<tor_err
        root=var_0;
        stepn=0;
    else
        [iter1,iter2]=nnl_newton_nle(nle_function,nle_derive,var_next,tor_err);
        root=1*iter1;
        stepn=1+iter2;
    end
end
