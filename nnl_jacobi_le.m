function solution_vec=nnl_jacobi_le(A_mat,b_vec,tor_err)

    f_norm2=@(vec)sqrt(sum(vec.^2));

    D_mat=diag(diag(A_mat));
    U_mat=triu(A_mat)-D_mat;
    L_mat=tril(A_mat)-D_mat;

    D_inv=inv(D_mat);

    B_jacobi=D_inv*(U_mat+L_mat);

    f_jacobi=D_inv*b_vec;

    x_k=repmat(0,[size(b_vec,1),1]);
    x_kp1=B_jacobi*x_k+f_jacobi;

    while(f_norm2(x_kp1-x_k)>tor_err)
        x_k=x_kp1;
        x_kp1=B_jacobi*x_k+f_jacobi;
    end
    
    solution_vec=x_kp1;
end

