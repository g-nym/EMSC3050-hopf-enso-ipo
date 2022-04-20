% Run to calculate a hopf bifurcation with noise.
 
close all
clear all
 
% Polar coordinates:
% dr/dt = lambda*r - mu*r^3
% dtheta/dt = omega
 
% Cartesian coordinates
% dx/dt = lambda*x - omega*y - mu*x*(x^2+y^2)
% dy/dt = lambda*y + omega*x - mu*y*(x^2+y^2)
 
% mu controls criticality of case (mu=1 - supercritical)



% write out noise vector
%dlmwrite('noise',randn(1,1e6),'-append','delimiter','\n');
 
% ENSO variables
% lambdae - self coupling strength/ damping
% omegae - frequency of ENSO
% sigmane - amplitude of noise added to ENSO
% sigmape - coupling strength of PDO to ENSO

% PDO variables
% lambdap - relaxation time
% omegap - frequency of PDO
% sigmanp - amplitude of noise added to PDO
% sigmaep - coupling strength of ENSO to PDO

lambdae = 0;
omegae = 2*pi/(4); 
noisee = 0.01; 
sigmape = 0.01; 
lambdap = 0; 
omegap = 2*pi/(4*15); 
noisep = 0.01; 
sigmaep = 0.01; 

if omegae==0 % checking if inherent timescale exists
    ite = 0;
    lambdae = 0;
else
    ite = 1;
end

if omegap==0
    itp = 0;
    lambdap = 0;
else
    itp = 1;
end

% Time variables
dt = 0.1;
stopt = 1000;
nt = stopt/dt;
 
% Initialise matrices
xe = zeros(2,nt);
dxe = zeros(1,2);
r2e = zeros(1,nt);
xp = zeros(2,nt);
dxp = zeros(1,2);
r2p = zeros(1,nt);

% Read noise
fid = fopen('noise');
dW = fscanf(fid,'%f',[2,nt-1]);
fclose(fid);
 
% Initial conditions
xe(1,1) = 0.01;
xe(2,1) = 0;
r2e(1) = (xe(1,1)*xe(1,1)) + (xe(2,1)*xe(2,1));

xp(1,1) = 0.01;
xp(2,1) = 0;
r2p(1) = (xp(1,1)*xp(1,1)) + (xp(2,1)*xp(2,1));

%  
for ii=2:nt
    % ENSO
    dxe(1) = lambdae*xe(1,ii-1) - omegae*xe(2,ii-1) - ite*xe(1,ii-1)*r2e(ii-1);
    dxe(2) = lambdae*xe(2,ii-1) + omegae*xe(1,ii-1) - ite*xe(2,ii-1)*r2e(ii-1);
    xe(:,ii) = xe(:,ii-1) + dxe(:)*dt + sigmape*xp(:,ii-1) + noisee*dW(:,ii-1);
    r2e(ii) = (xe(1,ii)*xe(1,ii)) + (xe(2,ii)*xe(2,ii)); 
    % PDO
    dxp(1) = lambdap*xp(1,ii-1) - omegap*xp(2,ii-1) - itp*xp(1,ii-1)*r2p(ii-1);
    dxp(2) = lambdap*xp(2,ii-1) + omegap*xp(1,ii-1) - itp*xp(2,ii-1)*r2p(ii-1);
    xp(:,ii) = xp(:,ii-1) + dxp(:)*dt + sigmaep*xe(:,ii-1) + noisep*dW(:,nt-ii+1);
    r2p(ii) = (xp(1,ii)*xp(1,ii)) + (xp(2,ii)*xp(2,ii));
end

taxis = 0:dt:(nt-1)*dt;

%%% Calculate spectra
yp = fft(xp(1,:));
yp(1) = [];
np = length(yp);
powerp = abs(yp(1:floor(np/2))).^2;
nyquist = 1/(dt*2); % half the sampling frequency
freqa = nyquist*linspace(0,1,(np/2));
periodp = 1./freqa; 
%perioda = 1./(2*pi*freqa)

ye = fft(xe(1,:));
ye(1) = [];
ne = length(ye);
powere = abs(ye(1:floor(ne/2))).^2;
nyquist = 1/(dt*2); % half the sampling frequency
freqe = nyquist*linspace(0,1,(ne/2));
periode = 1./freqe; 
%periode = 1./(2*pi*freqe)

%% Plot
 
tend = nt - 200/dt;

figure(1)
%afig(3)
subplot(211)
set(gca,'fontsize',12,'linewidth',1)
plot(taxis(tend:end),xe(1,tend:end),taxis(tend:end),xp(1,tend:end),'linewidth',1)
xlim([tend*dt stopt])
xlabel('time')
title('ENSO and PDO indices')
legend('ENSO','PDO','location','northeast')

subplot(223)
set(gca,'fontsize',12,'linewidth',1)
semilogx(periode,powere,'linewidth',1);
ylabel('Power');
xlabel('Period');
title('Spectrum of ENSO')
line([2*pi/omegae 2*pi/omegae], [0 1e6],'color','r')

subplot(224)
set(gca,'fontsize',12,'linewidth',1)
semilogx(periodp,powerp,'linewidth',1);
ylabel('Power');
xlabel('Period');
title('Spectrum of PDO')
line([2*pi/omegap 2*pi/omegap], [0 1e6],'color','r')


%print('-dpdf','simple.pdf')

