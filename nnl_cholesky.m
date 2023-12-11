function l_mat=nnl_cholesky(mat)
%%%CHOLESKY DECOMPOSITION FOR A QUARE MATRIX
%%% A=LL^T MAT=L_MAT*L_MAT^T

%% PRE MEMORIZE L MATRIX
    l_mat=zeros(size(mat));

    n=size(mat,1);

    %% PRE PROCCESS
    if  n==size(mat,2)
        if not(issymmetric(mat))
            error('ERROR:THE MATRIX MAT IS NOT SYMMETRIC!')
        else
            if not(all(eig(mat)>0))
                error('ERROR:THE MATRIX MAT IS NOT POSITIVE DEFINITE MATRIX!')
            end
        end
    else
        error('ERROR:THE MATRIX MAT MUST BE SQUARE!')
    end

    %% CHOLESKY DECOMPOSITION
    for j=1:n
        l_mat(j,j)=sqrt(mat(j,j)-sum(l_mat(j,1:j-1)).^2);
        for i=j+1:n
            l_mat(i,j)=(mat(i,j)-sum(l_mat(i,1:j-1).*l_mat(j,1:j-1)))/l_mat(j,j);
        end
    end
end
