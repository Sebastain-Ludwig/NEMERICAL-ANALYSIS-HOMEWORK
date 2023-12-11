function [q_mat,r_mat]=nnl_qr(mat)
%%%QR DECOMPOSITION FOR A QUARE(ILL CONDITION) MATRIX
%%% A=QR

%% 2-NORM
    f_norm2=@(vec)sqrt(sum(vec.^2));

    n=size(mat,1);
    if  n==size(mat,2)
        q_mat=q_de_recursive(n,mat);
    else
        error('ERROR:THE MATRIX MAT MUST BE SQUARE!')
    end

    
    %% PRE MEMORIZE Q R MATRIX
    function q_in=q_de_recursive(n,mat)
        if n==1
            q_in=[1];
        else
            omega=mat(:,1)-f_norm2(mat(:,1)).*eye(n,1);
            q_in=nnl_matexpand(q_de_recursive(n-1,f_matrix_ele(nnl_householder(omega)*mat,'2:end','2:end')),n)*nnl_householder(omega);
        end
    end

    r_mat=q_mat*mat;
end
