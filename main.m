%% PROBLEM 1-3:HORNER METHOD
ans_1_3=hornercal(ones(1,50),1.00001);
disp(["ANSWER 1-3:",ans_1_3]);

%% PROBLEM 2-1:LU DECOMPOSITION SOLVE A LINEAR EQUATION Ax=b
A_mat_2_1=[31,-13,0,0,0,-10,0,0,0;-13,35,-9,0,-11,0,0,0,0;0,-9,31,-10,0,0,0,0,0;0,0,-10,79,-30,0,0,0,-9;0,0,0,-30,57,-7,0,-5,0;0,0,0,0,-7,47,-30,0,0;0,0,0,0,0,-30,41,0,0;0,0,0,0,-5,0,0,27,-2;0,0,0,-9,0,0,0,-2,29];
lude(A_mat_2_1);

%% PROBLEM 2-2:CHOLESKY DECOMPOSITION SOLVE LINEAR EQUATION Ax=b
A_mat_2_2=[7,1,-5,1;1,9,2,7;-5,2,7,-1;1,7,-1,9];
b_vec_2_2=[13;-9;6;0];
l_mat_2_2=nnl_cholesky(A_mat_2_2);

%% PROBLEM 2-4:QR DECOMPOSITION
A_mat_2_4=[1 1 0 0;-1 3 -1/2 1/2;-2 2 3/2 1/2;-2 2 -1/2 5/2];
[q_mat_2_4,r_mat_2_4]=nnl_qr(A_mat_2_4);

%% PROBLEM 2-7:SOLVE LINEAR EQUATION (GAUSS_SEIDEL && JACOBI)
f_A_2_7_cons=@(n,m_var,l_var,u_var)diag(repmat(m_var,[1,n]))+diag(repmat(l_var,[1,n-1]),-1)+diag(repmat(u_var,[1,n-1]),1);

A_mat_2_7=arrayfun(f_A_2_7_cons,[10,20,30,50,100],repmat(3,[1,5]),repmat(-1,[1,5]),repmat(-1,[1,5]),'UniformOutput',false);
b_mat_2_7=arrayfun(@f_b_2_7_cons,[10,20,30,50,100],repmat(2,[1,5]),repmat(1,[1,5]),'Uniformoutput',false);

solution_vec_2_7=cellfun(@nnl_jacobi_le,A_mat_2_7,b_mat_2_7,{1e-5,1e-5,1e-5,1e-5,1e-5},'Uniformoutput',false);

%% PROBLEM 3-2:SOLVE NONLINEAR EQUATION
%% FUNCTION(A)
f_3_2_a=@(x)x*cos(x)+2;

root_3_2_a=nnl_dfr(-4,4,f_3_2_a,1e-3);

%% FUNCTION(B)
var_0_num=100;

f_3_2_b=cell(1,var_0_num);
fd_3_2_b=cell(1,var_0_num);

for i=1:var_0_num
    f_3_2_b{i}=@(x)x^3+2*x^2+10*x-100;
    fd_3_2_b{i}=@(x)3*x^2+2*x+10;
end

var_0_3_2_b=linspace(-10,20,var_0_num);
tor_err=repmat(1e-4,[1,var_0_num]);

[root_3_2_root,root_3_2_step]=cellfun(@nnl_newton_nle,f_3_2_b,fd_3_2_b,mat2cell(var_0_3_2_b,1,ones(1,var_0_num)),mat2cell(tor_err,1,ones(1,var_0_num)),'Uniformoutput',false);

step_3_2=[];
for i=1:var_0_num
    step_3_2=[step_3_2,root_3_2_step{i}(1)];
end
scatter(var_0_3_2_b,step_3_2);

%%PROBLEM 3-6:FIND ROOT ON A SPECIAL SPAN AND PROSIMATE ERROR
f_3_6=@(x)hornercal([-4,16,35,-69,-102,45,54],x);
fd_3_6=@(x)hornercal([16,35*2,-69*3,-102*4,45*5,54*6],x)

x_3_6=linspace(-2,2,100);
y_3_6=[];
y_3_6=arrayfun(f_3_6,x_3_6)

%%% PLOT
plot(x_3_6,y_3_6);

%%% FIND ROOTS
var_0_num=100;

cf_3_6=cell(1,var_0_num);
cfd_3_6=cell(1,var_0_num);

for i=1:var_0_num
    cf_3_6{i}=f_3_6;
    cfd_3_6{i}=fd_3_6;
end

var_0_3_6=linspace(-2,2,var_0_num);
tor_err=repmat(1e-4,[1,var_0_num]);

[root_3_6_root,root_3_6_step]=cellfun(@nnl_newton_nle,cf_3_6,cfd_3_6,mat2cell(var_0_3_6,1,ones(1,var_0_num)),mat2cell(tor_err,1,ones(1,var_0_num)),'Uniformoutput',false);

%%% GET UNIQUE ROOT FROM THE LIST LIST->SET
root_set=[]
for i=1:var_0_num
    root_set=[root_set,root_3_6_root{i}(1)];
end

root_set=uniquetol(root_set,1e-4);
