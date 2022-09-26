function A = afun_wrapper(i,j,x,dk,v0,dcorr)
% AFUN_HELM_SOUND_HARD_WRAPPER(I, J, X, ZPARS, NU, AREA, P, S) wraps the
% system matrix generator for the system of equations for the sound hard
% problem. I and J are index sets of the full set of points X for which we
% want to generate a system matrix. NU and AREA are the normal vectors and
% quadrature weights of the points X respectively. P is a permutation
% matrix required for skeletonization and S is a matrix of quadrature
% corrections.

% Convert system indices into original indices
ipts = idivide(int64(i(:)-1), int64(2))+1;
jpts = idivide(int64(j(:)-1), int64(2))+1;

[iuni,~,iiuni] = unique(ipts);
[juni,~,ijuni] = unique(jpts);

% Get matrix valued entries of system matrix for unique points
A_uni = afun(iuni,juni,x,dk,v0,dcorr);

% Extract relevant rows and columns in A_uni
iiuni2 = (iiuni-1)*2 + mod(i(:)-1, 2)+1;
ijuni2 = (ijuni-1)*2 + mod(j(:)-1, 2)+1;

A = A_uni(iiuni2, ijuni2);
end

