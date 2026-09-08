%Combine the basis for the tangent space of the T_N manifold with the basis for the tangent space of \mathcal{K}_{F,\eta,q}

%Run this script only after symbolically computing the T_N configuration,
%resulting in the X symbolic variable being in the MATLAB workspace.

for I=1:N

     %we cyclically shift the indices in order to check that T_X K +
    %ker D\pi_k=(\mathbb{R}^{3\times2})^N for each k=1,...,N
  %FIRST, define a new starting point for the T_N configuration
P=P+C(:,:,1);

%now, do the circular shifting
C=circshift(C,-1,3);
kappa=circshift(kappa,-1,1);
X=circshift(X,-1,3);



%first, need to divide the rank-one matrices C(:,:,i) into 2-dimensional
%and 3-dimensional vectors

%first, factorize the rank-one matrices C(:,:,i) into their constitutive
%2-dimesional and 3-dimensional vectors
[aCell, bCell] = factorRankOne3x2Symbolic(C);

out = tn_tangent_basis_symbolic(P, aCell, bCell, kappa);
%generate basis for tangent space of T_N manifold
%deleting the first 6 vectors gives a basis for the kernel of the derivative of \pi
kernel_Pi_basis=out.basisMatrix(:,7:end);

%now, calculate the basis for the tangent space of \mathcal{K}_{f,\eta,q}
basisK=[];
basis_state_space %script that populates basisK with the tangent space of \mathcal{K}_{F,\eta,q}

%candidate for the final basis
basis=[kernel_Pi_basis basisK];
disp('rank of potential basis set is')
rank(basis)
end

%note that after running the above loop, the $T_N$ is back in the original
%order
