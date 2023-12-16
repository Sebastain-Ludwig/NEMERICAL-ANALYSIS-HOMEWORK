function [l_mat,u_mat]=lude(mat)
%%%DOOLITTLE DECOMPOSITION FOR A QUARE MATRIX
%%% A=LU MAT=L_MAT*U_MAT

    
%% PRE MEMORIZE L U MATRIX
	u_mat=zeros(size(mat));
	l_mat=zeros(size(mat));

        %% PRE CALCULATE L U FISRT COL %% ROW
	u_mat(1,:)=mat(1,:);

	zerop_f=@(mat_ele) (mat_ele>1e-4);

	l_mat(:,1)=mat(:,1)./mat(1,1);

	n=size(mat,1);

        %% DOLITTILE MATHED CALCULATE
	for i=2:n
		for j=i:n
			u_mat(i,j)=mat(i,j)-l_mat(i,1:i-1)*u_mat(1:i-1,j);
			if zerop_f(u_mat(i,i))
				l_mat(j,i)=(mat(j,i)-l_mat(j,1:i-1)*u_mat(1:i-1,i))/u_mat(i,i);
			else
				error('ERROR:ELEMENR MATTIRX MAIN NUMERICAL LOSS,LOCATION:%s,%s',num2str(i),num2str(j));
			end
		end
	end
end
