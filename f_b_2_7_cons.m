function cell_b=f_b_2_7_cons(n,side_var,mid_var)
    mat_end_ele=zeros(n,1);
    mat_end_ele(end)=side_var-mid_var;
    cell_b=repmat(mid_var,[n,1])+(side_var-mid_var)*eye(n,1)+mat_end_ele;
end
