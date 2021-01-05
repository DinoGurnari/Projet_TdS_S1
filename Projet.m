F0 = 6000; %Hz
F1 = 2000; %Hz
Fe = 48000; %Hz

Te = 1/Fe; %s
Ts = 1/300; %s/bits
Ns = floor(Ts/Te); %échantillons/bits

donnees = randi(0:1, 30, 1);
NRZ = zeros(length(donnees)*Ns,1);
T = ([0:Ns*length(donnees)-1]*Te)';

%3.1.1.1

for i = 1:length(donnees)
    NRZ((i-1)*Ns+1:i*Ns) = donnees(i);
end;

%3.1.1.2

figure;
plot(T,NRZ);
xlabel('t en s');
ylabel('NRZ(t)');

%3.1.1.3

S_NRZ = abs(fft(NRZ).^2)/length(NRZ);
figure;
subplot(211);
semilogy(S_NRZ);
subplot(212);
pwelch(NRZ);

%3.1.2.1

phi0 = rand*2*pi;
phi1 = rand*2*pi;

Cos0 = cos(2*pi*F0*T + phi0);
Cos1 = cos(2*pi*F1*T + phi1);
Un = eye(size(NRZ));

x = (1 - NRZ).*Cos0 + NRZ.*Cos1;

%3.1.2.2

figure;
plot(T,x);

%3.1.2.3

%3.1.2.4

S_x = abs(fft(x).^2)/length(x);
figure;
subplot(411);
semilogy(S_x);
subplot(412);
pwelch(x);

%3.3.1 Filtre passe-bas

g = Fe*(1:length(x))/length(x);
fc = 2500; %Hz
intervalle = [-100*Te:Te:100*Te];

figure;
h = (2*fc/Fe)*sinc(2*fc*intervalle);
H = 50*fft(h);
z = filter(H,[1],x);
Z = abs(fft(z));
plot(intervalle,H);

%3.3.2 Filtre passe-haut
