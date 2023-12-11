function matrix_ctr=f_matrix_ele(matrix_in,row_str,col_str)
    mt_in=matrix_in;
    matrix_ctr=eval(['mt_in','(',row_str,',',col_str,')']);
end
