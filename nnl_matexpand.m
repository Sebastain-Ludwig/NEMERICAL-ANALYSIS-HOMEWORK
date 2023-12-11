function expand_mat=nnl_matexpand(mat,exp_size)
    n=size(mat,1);
    if  n==size(mat,2)
        if exp_size<n
            error('ERROR:THE EXPAND MATRIX SIZE IS SMALLER THAN THE MAT SIZE!')
        end
    else
        error('ERROR:THE MATRIX MAT MUST BE SQUARE!')
    end

    expand_mat=zeros(exp_size);
    expand_mat(exp_size-n+1:end,exp_size-n+1:end)=mat;
    expand_mat(1:exp_size-n,1:exp_size-n)=eye(exp_size-n);
end
