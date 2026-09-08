%This script calculates the basis for the tangent space of \mathcal{K}_{f,\eta,q}
for i=1:2*N
    element=sym(zeros(2*N,1));
    element(i)=1;
    basis_element=sym(zeros(3,2,N));
    DF=[0 -1;-1 0]; %choose p'=-1
    for k=1:N
        u=element((2*(k-1)+1):(2*(k-1)+2),1); %the parameter we use
        basis_element(1:2,1,k)=u;
        basis_element(1:2,2,k)=-DF*u;
        basis_element(3,1,k)=[X(2,2,k) X(2,1,k)]*u;
        basis_element(3,2,k)=-[X(2,1,k)*(-1) -X(2,2,k)]*u; %choose p'=-1
    end
        basisK(:,end+1)=reshape(basis_element,[6*N,1]); %turn each basis element into one big vector
end