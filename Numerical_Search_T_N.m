%This is the MATLAB code to numerically find a potential T_N configuration
%for the p-system, in the space of 3x2 matrices

exitflag=0;

%keep looping until the solver has found a solution
while exitflag~=1

%clear the workspace variables (this is necessary)
clear all

%the epsilon determines the strictness in strict inequalities
epsilon=10;

%pick N\geq4 for the T_N configuration we are looking for
N=6;
%kappa=optimvar('kappa',N,1,LowerBound=epsilon);
n=optimvar('n',1,2,N);
a=optimvar('a',3,1,N);
X=optimvar('X',3,2,N);
P=optimvar('P',3,2);

%the kappa_i are fixed
kappa=rand(N,1)*1000+1;


prob=optimproblem;


C=optimexpr(3,2,N);
for i=1:N
C(:,:,i)=a(:,:,i)*n(:,:,i);
end

%we want the C_i to sum to zero (so sum along the third component of the C
%array)

prob.Constraints.C_sum_to_zero=sum(C,3)==zeros(3,2);


%we want a convex entropy (i.e., this forces p'<0)
for i=1:N
for j=[1:(i-1) (i+1):N]
constraint = strcat('ConvexityCondition','_',num2str(i),'_',num2str(j));
prob.Constraints.(constraint)=(X(3,1,j)-.5*X(2,1,j)^2)-(X(3,1,i)-.5*X(2,1,i)^2)-X(2,2,i)*(X(1,1,j)-X(1,1,i))>=epsilon;
end
end

for i=1:N
%we want the following symmetry-type condition on the T_N
constraint = strcat('SymmetryCondition1','_',num2str(i));
prob.Constraints.(constraint)=X(2,1,i)-X(1,2,i)==0;

%force the condition on q

constraint = strcat('qCondition','_',num2str(i));
prob.Constraints.(constraint)=X(3,2,i)-X(1,2,i)*X(2,2,i)==0;
end



% we want the X to be in T_N configuration
for i=1:N
constraint = strcat('constrT_N_Configuration',num2str(i));
prob.Constraints.(constraint)=X(:,:,i)-P-sum(C(:,:,1:(i-1)),3)-kappa(i)*C(:,:,i)==zeros(3,2);
end

%initialize the solver with random initial points
data.a = randi([-1000,1000],3,1,N);
data.n=randi([-1000,1000],1,2,N);
data.X=randi([-1000,1000],3,2,N);
data.P=randi([-1000,1000],3,2);
%data.kappa=randi([1,100],N,1);
%or, use previously found (symbolic) values as a starting point:
%data.a = double(sol.a);
%data.n=double(sol.n);
%data.X=double(sol.X);
%data.P=double(sol.P);


%solve!
disp('Start the solver')
%[sol,fval,exitflag,output] = solve(prob,data,'Options',optimoptions(prob,'MaxFunctionEvaluations',200000,'MaxIterations',20000000,'ConstraintTolerance',.00000000000001,'EnableFeasibilityMode',true,'SubproblemAlgorithm','cg','Display','iter'));
[sol,fval,exitflag,output] = solve(prob,data,'Options',optimoptions(prob,'MaxFunctionEvaluations',200000,'MaxIterations',50000,'ConstraintTolerance',.00000000000001,'Display','iter'));

end

%convert numeric values to symbolic expressions
sol.a=sym(sol.a,'f');
sol.n=sym(sol.n,'f');
sol.P=sym(sol.P,'f');
%kappa=sym(sol.kappa,'f');
kappa=sym(kappa,'f');
%create the sum of the C_i symbolically
S=sym(zeros(3,2,1));
for i=1:N
S=S+sol.a(:,:,i)*sol.n(:,:,i);
end
