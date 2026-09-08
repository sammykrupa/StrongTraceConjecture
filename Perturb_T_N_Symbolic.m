%This is the script to perturb a numerical solution to the T_N equations
%found using Numerical_Search_T_N.m into an exact solution which
%(symbolically) verifies all of the necessary conditions

%This code assumes that the script Numerical_Search_T_N.m has succeeded in finding a numerical solution for a T_N
% with N=6, or the data has been loaded from a MAT file back into MATLAB, and the MATLAB workspace contains the relevant variables
N=6;

%save original sol structure from the gradient descent search
solOriginal=sol;

%calculate the Xs and Cs
X=sym(zeros(3,2,N));
C=sym(zeros(3,2,N));
for i=1:N
C(:,:,i)=sol.a(:,:,i)*sol.n(:,:,i);
end
S=sym(zeros(3,2));
for i=1:N
X(:,:,i)=sol.P+S+kappa(i)*C(:,:,i);
S=S+C(:,:,i);
end

%FORCE the symmetry-type condition to be EXACTLY true
for i=1:N
X(2,1,i)=X(1,2,i);
end
solOriginal.P(2,1)=solOriginal.P(1,2);

%save the original X values
Xold=X;

X=X(1:2,1:2,1:N);

%re-using code from the paper "Non-Uniqueness and BV Blowup for Vacuously Liu-Admissible Solutions to P-System Via Computer-Assisted Proof" by Krupa (J. Math. Fluid Mech., 28(2):Paper No. 37, 2026)
%This code assumes that different orderings of the matrices in the T_N
%configuration each lead to a different parameterization of the points in the
%T_N. However, here we restrict the code to the single ordering [1:N].
%we compute the paramtrization (P,C_i,kappa_i) (following
%Székelyhidi, László, Jr.: Rank-one convex hulls in R^{2×2}. Calc.Var.PDE 22(3),253–281(2005))
%see also Förster, Clemens and Székelyhidi, László, Jr.: T5-configurations and non-rigid sets of matrices. Calc. Var. Partial Differential Equations 57(1), Art. 19 (2018)
%We use the ordering 1:N
SavedOrderings=[1:N]; %one ordering for each row, but we only use this one ordering.
%initialize the symbolic matrices we will use to store our
%parameterizations
sol.P=sym(zeros(2,2,3));
sol.C=sym(zeros(2,2,N,3));
sol.kappa=sym(zeros(N,1,3));

for q=1:1 %loop through the orderings, but we only consider 1

%symbolically find the scalar mu and the vector lambda 

syms mu;

    for i=1:N
    for j=1:N
        temp=X(:,:,SavedOrderings(q,i))-X(:,:,SavedOrderings(q,j));
        if i>j
            %compute determinant of temp matrix (multiplied by mu)
            Amu(i,j)=mu*(det(temp));
        else
            %compute determinant of temp matrix
            Amu(i,j)=det(temp);
        end
    end
    end

%solve for mu symbolically using root function 
%see https://www.mathworks.com/help/symbolic/sym.root.html

zerosVal=root(det(Amu),mu);
muValue=zerosVal(5); %for the given orderings in SavedOrdering, the fifth zero always give an acceptable value of mu (>1), at least for the T_6 example from the paper itself
%substitute this value of mu into the matrix Amu, then calculate a basis
%for the null space of Amu
%basis=null(subs(Amu,mu,muValue)); %for the given orderings in SavedOrdering, this basis is one vector, each component of which is positive. Thus, we can take basis=lambda;
%However, the null command does not always work because it needs to simplify
%complicated symbolic expressions. Instead, we manually find the
%eigenvector by using that each nonzero column of the adjugate matrix is in
%the null space. From ChatGPT:
A = Amu;

AdjA = sym(zeros(N,N));
for i = 1:N
    for j = 1:N
        M = A;
        M(i,:) = [];
        M(:,j) = [];
        AdjA(j,i) = (-1)^(i+j) * det(M); 
        % note transpose: adj(A) is transpose of cofactor matrix
    end
end

AdjAtRoot = subs(AdjA, mu, muValue);

lambda = sym([]);
for j = 1:N
    col = AdjAtRoot(:,j);
    if ~isAlways(all(col == 0))
        lambda = col;
        break
    end
end

if isempty(lambda)
    error('All columns of adj(A(muValue)) vanished; rank may be < N-1 or the chosen root is wrong.');
end
 

%Compute the parameterization of the T_N. 

%Now, from the algebraic criterion, we can compute the parameterization
%(P,C_i,kappa_i) for this ordering of the T_N

 C=sym(zeros(2,2,N));
 kappa=sym(zeros(N,1));
 tempP=sym(zeros(2,2,N)); %for each ordering of the T_N, this will store the P_i associated with the vertices of the inner `rank-one' loop


    %compute the P_i
    for i=1:N
            %compute the scalar pre-factor
            denominator=sum(lambda(i:N,1))+muValue*sum(lambda(1:(i-1),1));
            %compute the matrix term
            matrixSum=sym(zeros(2,2));
            for k=1:N
                if k<i
                    matrixSum=matrixSum+muValue*lambda(k,1)*X(:,:,SavedOrderings(q,k));
                else
                    matrixSum=matrixSum+lambda(k,1)*X(:,:,SavedOrderings(q,k));
                end
            end
   tempP(:,:,i)=(1/denominator)*matrixSum;

         if i==1
           sol.P(:,:,q)=tempP(:,:,i); %save the initial P values for the first, second and third orderings. However, here we are only concerned with the first ordering
         end
         
         if i>1
             tempC=tempP(:,:,i)-tempP(:,:,i-1); %this gives the C_{i-1}
             tempkappaC=X(:,:,SavedOrderings(q,i-1))-tempP(:,:,i-1); %this gives kappa_i C_{i-1}
             kappa(i-1,1)=tempkappaC(1,1)/tempC(1,1); %assuming tempC(1,1) is not zero
             C(:,:,i-1)=tempC; %save the C_i
         end

    end
    %now, collect the info from the last of the P_i (i=N) after the for
    %loop above has finished
             tempC=tempP(:,:,1)-tempP(:,:,N); %this gives the C_N
             tempkappaC=X(:,:,SavedOrderings(q,N))-tempP(:,:,N); %this gives kappa_N C_N
             kappa(N,1)=tempkappaC(1,1)/tempC(1,1);
             C(:,:,N)=tempC; %save the C_N

              %save all of our results from determining the
              %parameterization corresponding to this particular ordering
              %of the T_N (ordering #q)
             sol.C(:,:,:,q)=C;
             sol.kappa(:,:,q)=kappa;
end

%now, work on third row of the T_N
%want to ensure the third row of the C_i sums to zero


 %calculate new n vectors
for i=1:N
   sol.n(1,2,i)=sol.n(1,1,i)*(C(1,2,i)/C(1,1,i));
end

%build third row
for i=1:N
   C(3,:,i)=sol.n(1,:,i)*sol.a(3,1,i);
end

%calculate sum of first N-2 rank-one directions
partial_loop=-(sum(C(3,:,1:(N-2)),3))';

%some linear algebra
n_matrix=[sol.n(1,1,5) sol.n(1,1,6);sol.n(1,2,5) sol.n(1,2,6)];
a_matrix=inv(n_matrix)*partial_loop;

%calculate new sol.a values
sol.a(3,1,N-1)=a_matrix(1);
sol.a(3,1,N)=a_matrix(2);

%rebuild third row
for i=1:N
   C(3,:,i)=sol.n(1,:,i)*sol.a(3,1,i);
end

%calculate the Xs
X=sym(zeros(3,2,N));

S=zeros(3,2);
for i=1:N
X(:,:,i)=solOriginal.P+S+kappa(i)*C(:,:,i);
S=S+C(:,:,i);
end

%now, calculate the perturbations of the kappas in order to get the correct
%values in the (3,2) slot of the T_N

syms kappa_var [1 6]
S=sym(zeros(3,2));
for i=1:N
X(:,:,i)=solOriginal.P+S+kappa_var(i)*C(:,:,i);
S=S+C(:,:,i);
end

for i=1:N
    i
    kappa_sol=solve(X(3,2,i)-X(1,2,i)*X(2,2,i));
    [~,I]=min([double(abs(kappa_sol(1)-kappa(i))),double(abs(kappa_sol(2)-kappa(i)))]);

    if I==1
        kappa(i)=kappa_sol(1);
    else
        kappa(i)=kappa_sol(2);
    end
end


%do things by hand
% for i=1:N
%     i
%     disp('original error')
%     error=(X(3,2,i)-X(1,2,i)*X(2,2,i));
%     vpa(error)
%     AKold=kappa(i)*C(3,2,i)-kappa(i)^2*C(2,2,i)*C(1,2,i);
% 
%     %compute discriminant
%     disc = C(3,2,i)^2 - 4*(-C(2,2,i)*C(1,2,i))*(error - AKold);
% 
%     KappaNewPlus  = (-C(3,2,i) + sqrt(disc)) / (2*(-C(2,2,i)*C(1,2,i)));
%     KappaNewMinus = (-C(3,2,i) - sqrt(disc)) / (2*(-C(2,2,i)*C(1,2,i)));
% 
%     %vpa(kappa(i))
%     %vpa(KappaNewPlus)
%     %vpa(KappaNewMinus)
% 
%     [~,I]=min([double(abs(KappaNewPlus-kappa(i))),double(abs(KappaNewMinus-kappa(i)))]);
%     I
%     if I==1
%         kappa(i)=KappaNewPlus;
%     else
%         kappa(i)=KappaNewMinus;
%     end
% end

%now, rebuild the T_N
X=sym(zeros(3,2,N));

S=sym(zeros(3,2));
for i=1:N
X(:,:,i)=solOriginal.P+S+kappa(i)*C(:,:,i);
S=S+C(:,:,i);
end




%%%SYMBOLIC VERIFICATIONS OF ALL REQUIREMENTS%%%%
%We have a rank-one loop
disp('compute rank of C matrices')
for i=1:N
    A=C(:,:,i);
   %we do this by hand, because rank(A) sometimes does not like how complicated the symbolic expressions are becoming
 
% Rank-one test:
% all 2x2 minors must be zero, and A must not be the zero matrix.
%Sometimes MATLAB throws a warning, because things are too complicated for
%the simplify command to workout. However, this rank-one test is not really
%necessary. The Székelyhidi theory always ensures that we have a rank-one loop.
m1 = simplify(det(A([1 2], :)))
m2 = simplify(det(A([1 3], :)))
m3 = simplify(det(A([2 3], :)))

isRankOne = isAlways(m1 == 0) && ...
            isAlways(m2 == 0) && ...
            isAlways(m3 == 0) && ...
            ~all(arrayfun(@(x) isAlways(x == 0), A(:)));

disp(isRankOne)
end

%confirm we really have a loop
isAlways(simplify(sum(C,3))==sym(zeros(3,2)))

%kappa values always greater than one
disp('kappa values need to be greater than one')
vpa(kappa)

%condition on the entropy-flux 
disp('Check condition on the entropy-flux')

for i=1:N
    isAlways(simplify(X(3,2,i)-X(1,2,i)*X(2,2,i))==0)
end

%symmetry condition
disp('Check symmetry condition')

for i=1:N
    isAlways(X(2,1,i)==X(1,2,i))
end

%strict convexity of the entropy
disp('Checking convexity of entropy -- want these greater than 0')
for i=1:N
for j=[1:(i-1) (i+1):N]
    vpa((X(3,1,j)-.5*X(2,1,j)^2)-(X(3,1,i)-.5*X(2,1,i)^2)-X(2,2,i)*(X(1,1,j)-X(1,1,i)))
end
end