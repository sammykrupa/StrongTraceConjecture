%This script checks symbolically that the T_6 from ChatGPT meets our
%requirements.
N=6;
data=T6_full_data(true);
C=sym(zeros(3,2,6));
for i=1:N
    C(:,:,i)=sym(data.C(i));
end
kappa=sym(data.kappa);
P=sym(data.P);

%make the X matrices
X=sym(zeros(3,2,N));

S=sym(zeros(3,2));
for i=1:N
X(:,:,i)=P+S+kappa(i)*C(:,:,i);
S=S+C(:,:,i);
end




%%%SYMBOLIC VERIFICATIONS OF ALL REQUIREMENTS%%%%
%We have a rank-one loop
disp('compute rank of C matrices')
for i=1:N
    rank(C(:,:,i))
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

%check no rank-one connections
disp('Dont want rank-one connections between the T6 elements, so these ranks should all be 2')

for i=1:N
for j=[1:(i-1) (i+1):N]
rank(X(:,:,i)-X(:,:,j))
end
end

