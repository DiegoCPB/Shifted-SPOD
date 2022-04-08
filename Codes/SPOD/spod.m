function [Psi, Lambda, Qhat, St, Nb] = spod(Q,W,dt,nfft,olap)
% Dimensions
% ----------
% N :                     number of spatial points, d.o.f.
% Nt :                    number of snapshots/realizations. 
% NSt :                   number of frequencies.
% Nb :                    number of Welch blocks.
% Nm :                    number of SPOD modes. min(N,Nb).

% INPUTs
% ------
% Q [N x Nt] :            Data matrix, without mean.
% W [N x 1] or [1 x N] :  Integration weights.
% dt [1 x 1] :            Time step.
% nfft [1 x 1]:           Block size for Welch method.
% olap [1 x 1]:           Overlap in % for Welch method.
% Default window is infinitely continuous (https://arxiv.org/pdf/1907.04787.pdf).

% OUTPUTs
% -------
% Psi [NSt x N x Nm] :    Spatial SPOD modes
% Lambda [NSt x Nm] :     SPOD gains
% Qhat [NSt x N x Nb] :   Estimated spectral data matrices
% St [1 x NSt] :          Non-dimensional frequency

%========================================================
% Authors : Leandra I. Abreu, Diego C. P. Blanco 
% Last version : 08/04/2022
% Referencias : 
% Abreu, 2017
% Cordier e Bergmann, 2003 (POD) 
% Towne 2018, (POD, SPOD)
% Schimdt 2020
%========================================================

[N, Nt] = size(Q);
Nb = calcNb(Nt,nfft,olap);
fs = 1/dt; % sampling frequency

St = (fs/nfft)*(0:(nfft-1));
St(St >= fs/2) = St(St >= fs/2) - fs; % two sided spectrum
NSt = length(St);   

Qhat = zeros(NSt,N,Nb);

% Windowing
w = inf_smooth(nfft);
wrms = rms(w);
ECF = 1/wrms; % Window energy correction factor
Ww = repmat(w,[N 1]);

for i = 1:Nb
    pos1 = calcPos1(olap,nfft,i);   
    pos2 = pos1+nfft-1;
    Qfft = ECF*fft(Q(:,pos1:pos2).*Ww,[],2)/nfft;
    Qhat(:,:,i) = Qfft.';
end

Nm = min(N,Nb);
Psi = zeros(NSt,N,Nm);
Lambda = zeros(NSt,Nm);

W = spdiags(W(:),0,N,N);
for j=1:NSt
    Qhat_aux = squeeze(Qhat(j,:,:));
    
    % Direct POD
    if Nm == N
        C = Qhat_aux*Qhat_aux'*(W/Nb); 
        [Psi_aux,lambda] = eig(C);
        [lambda,ind_sort] = sort(diag(lambda),'descend');
        Psi_aux = Psi_aux(:,ind_sort); 
        Psi_aux = Psi_aux./sqrt(sum(Psi_aux.*conj(Psi_aux).*diag(W),1)); % Normalisation
        
    % Snapshot POD
    elseif Nm == Nb
        C = Qhat_aux'*(W/Nb)*Qhat_aux;
        [Theta_aux,lambda] = eig(C);
        [lambda,ind_sort] = sort(diag(lambda),'descend');
        Theta_aux = Theta_aux(:,ind_sort); 
        Psi_aux = Qhat_aux*Theta_aux./sqrt(Nb*lambda.'); 
    end

    Psi(j,:,:) = Psi_aux; % spatial modes
    Lambda(j,:) = lambda.';
end

end

function Nb = calcNb(Nt,Nfft,olap)
    Nolap = floor(Nfft*olap/100);
    Nb = floor((Nt-Nolap)/(Nfft-Nolap)); % [Schmidt, 2020], pag 6
end

function pos1 = calcPos1(olap,nfft,i)
    pc = olap/100;
    pos1 = (i-1)*floor(nfft*(1-pc))+1;
end

function y = inf_smooth(n)
    x = linspace(0,1,n);
    y = zeros(1,n);
    y(2:end-1) = exp(8)./exp(2./(x(2:end-1).*(1-x(2:end-1))));
end