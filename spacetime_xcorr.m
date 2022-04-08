clear all

solution = matfile('GL_solution.mat');
q = solution.q;
x = solution.x;
t = solution.t;

U = 10;
xp = 15;
nfft = 100;
n = size(q,1);
dt = t(2)-t(1);

[~,idx] = min(abs(xp-x));
xp = x(idx);

shift = (xp-x)./(U*dt);
qps = q;
for i = 1:size(qps,1)
    S = shift(i);
    Sfloor = floor(S);
    Sceil = ceil(S);
    if S ~= Sfloor  % Linearly interpolated shift
        qps(i,:) = (Sceil-S)*circshift(qps(i,:),Sfloor)+(S-Sfloor)*circshift(qps(i,:),Sceil);
    else % Direct shift
        qps(i,:) = circshift(qps(i,:),Sfloor);
    end
end

ref_signal = q(idx,:);
cross_corr = [];
for ii=1:n
    signal = q(ii,:);
    [cross, lags] = xcorr(ref_signal,signal,max(abs(shift)));
    cross_corr = [cross_corr; cross];
end

ref_signal = q(idx,:);
cross_corr2 = [];
for ii=1:n
    signal = qps(ii,:);
    [cross, lags] = xcorr(ref_signal,signal,max(abs(shift)));
    cross_corr2 = [cross_corr2; cross];
end

[X,L] = meshgrid(xp-x,lags*dt);
C = transpose(abs(cross_corr));
C2 = transpose(abs(cross_corr2));

set(0,'defaultTextInterpreter','latex'); 

figure();
subplot(1,2,1)
hold on;
contourf(L,X,C,'LineColor','none');
xline(-nfft*dt/2); xline(nfft*dt/2);
ylabel('$\Delta x$')
xlabel('$\Delta t$')
title('Space-time correlations')

subplot(1,2,2)
hold on;
contourf(L,X,C2,'LineColor','none');
colorbar();
xline(-nfft*dt/2); xline(nfft*dt/2);
xlabel('$\Delta t$')
title('Shifted correlations')