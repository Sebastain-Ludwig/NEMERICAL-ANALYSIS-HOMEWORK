function hs_mat=nnl_householder(omega)
    if size(omega,2)~=1
        error('OMEGA VEC MUST BE COL VEC!')
    end
    hs_mat=eye(size(omega,1))-(2/(transpose(omega)*omega)).*(omega*transpose(omega));
end
