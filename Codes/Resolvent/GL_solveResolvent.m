function [U,S,V,x] = GL_solveResolvent(LinearMatrixPath,omega)
    LINEAR = load(LinearMatrixPath);
    L = LINEAR.L;
    x = LINEAR.x;

    II = eye(size(L));

    % Weights
    W = weightVector(x);
    W = diag(W);

    % Resolvent operator
    R = -inv(1i*omega*II - L);

    % optimal forcing according to Rayleigh quotient:
    % max <Rf,Rf>/<f,f>
    % max fH RH W1/2H W1/2 R f / (fH W1/2H W1/2 f)
    % max gH W-1/2H RH W1/2H W1/2 R W-1/2 g / (gH g) with g=W1/2 f

    Wchol = sqrtm(W);

    [U,S,V]=svd(Wchol*R/Wchol);

    %see Lesshafft, Semeraro, Jaunet, Cavalieri & Jordan (arXiv 2018):
    %Wchol*R*inv(Wchol) =U Sigma V';
    %R =inv(Wchol)*U Sigma V'*Wchol;
    %optimal forcing/response given by inv(Wchol)*V and inv(Wchol)*U

    V = Wchol\V;
    U = Wchol\U;
end