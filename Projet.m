clear all
close all
clc

% Paramètres
F0 = 6000; % Fréquence bits 0 en Hz
F1 = 2000; % Fréquence bits 1 en Hz
Fe = 48000; % Fréquence d'échantillonage en Hz
Te = 1/Fe; % Période d'échantillonage en secondes
Ts = 1/300; % s/bits, le débit souhaité max est de 300 bits/seconde

Ns = floor(Ts/Te); % échantillons/bits, on prend la partie entière, Ts = Ns*Te
Nb_bits = 200; % Nombre de bits du signal émis
Nb_echantillons = Nb_bits*Ns; % Nombre d'échantillons

donnees = randi(0:1, Nb_bits, 1); % Génération de Nb_bits 0 ou 1 de manière aléatoire
NRZ = zeros(Nb_echantillons,1); % Initialisation du signal NRZ avec Nb_echantillons echantillons
T = ([0:Nb_echantillons-1]*Te)'; % Le signal est tracé de 0 à (Nb_echantillons-1)*Te secondes,
                                 % avec un pas de Te

                                 
%% Modem de fréquence - Démodulation par filtrage
        
%3.1 Construction du signal modulé en fréquence

%3.1.1.1 Génération du signal NRZ 
for i = 1:Nb_bits
    NRZ((i-1)*Ns+1:i*Ns) = donnees(i); 
end;

%3.1.1.2 Tracé du signal NRZ(t) en fonction du temps
figure;
plot(T,NRZ);
ylim([-0.5 1.5]);
xlabel('t en s');
ylabel('NRZ(t)');
title('Signal NRZ(t) en fonction du temps');

%3.1.1.3 Tracé de la DSP de NRZ

[DSP_NRZ_Th,f] = DSP_rectangulaire(Nb_echantillons,Fe,NRZ);
figure;
subplot(211);
semilogy(f,DSP_NRZ_Th);
title('DSP théorique de NRZ en fonction de f');

DSP_NRZ_Exp = (Ts*sinc(pi*f*Ts).^2 + dirac(f))/4;
%figure;
subplot(212);
semilogy(f,DSP_NRZ_Exp);
title('DSP expérimentale de NRZ en fonction de f');

%DSP_signal = periodogram(signal,[],Nb_points,'centered');
figure;
pwelch(NRZ);

% %3.1.2.1 Génération du signal x(t) modulé en fréquence
% 
% phi0 = rand*2*pi;
% phi1 = rand*2*pi;
% 
% Cos0 = cos(2*pi*F0*T + phi0);
% Cos1 = cos(2*pi*F1*T + phi1);
% Un = eye(size(NRZ));
% 
% x = (1 - NRZ).*Cos0 + NRZ.*Cos1;
% 
% %3.1.2.2
% 
% figure;
% plot(T,x);
% 
% %3.1.2.3
% 
% %3.1.2.4
% 
% S_x = abs(fft(x).^2)/length(x);
% figure;
% subplot(411);
% semilogy(S_x);
% subplot(412);
% pwelch(x);
% 
% %3.3.1 Filtre passe-bas
% 
% g = Fe*(1:length(x))/length(x);
% fc = 2500; %Hz
% intervalle = [-100*Te:Te:100*Te];
% 
% figure;
% h = (2*fc/Fe)*sinc(2*fc*intervalle);
% H = 50*fft(h);
% z = filter(H,[1],x);
% Z = abs(fft(z));
% plot(intervalle,H);
% 
% %3.3.2 Filtre passe-haut
