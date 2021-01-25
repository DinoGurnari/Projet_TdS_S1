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
Nb_bits = 100; % Nombre de bits du signal émis
Nb_echantillons = Nb_bits*Ns; % Nombre d'échantillons
Null = 2000;

donnees = randi(0:1, Nb_bits, 1); % Génération de Nb_bits 0 ou 1 de manière aléatoire
NRZ = zeros(Nb_echantillons,1); % Initialisation du signal NRZ avec Nb_echantillons echantillons
T = ([0:Nb_echantillons-1]*Te)'; % Le signal est tracé de 0 à (Nb_echantillons-1)*Te secondes,
                                 % avec un pas de Te

                                 
%% 3.Modem de fréquence - Démodulation par filtrage
        

%3.1 Construction du signal modulé en fréquence

%3.1.1.1 Génération du signal NRZ 
for i = 1:Nb_bits
    NRZ((i-1)*Ns+1:i*Ns) = donnees(i); 
end;

%3.1.1.2 Tracé du signal NRZ(t) en fonction du temps
figure; % figure 1
plot(T,NRZ);
ylim([-0.5 1.5]);
xlabel('t en s');
ylabel('NRZ(t)');
title('Signal NRZ(t) en fonction du temps');

%3.1.1.3 Tracé de la DSP de NRZ théorique et expérimentale
[DSP_NRZ_Th,f] = DSP_rectangulaire(Nb_echantillons,Fe,NRZ,0);
DSP_NRZ_Exp = (Ts*sinc(pi*f*Ts).^2 + dirac(f))/4;

figure; % figure 2
subplot(211);
semilogy(f,DSP_NRZ_Th);
ylim([1e-06 1000]);
xlim([-2000 2000]);
xlabel('f en Hz');
ylabel('DSP NRZ Th(f)');
title('DSP théorique de NRZ en fonction de f');

subplot(212);
semilogy(f,DSP_NRZ_Exp);
ylim([1e-10 1e-2]);
xlim([-2000 2000]);
xlabel('f en Hz');
ylabel('DSP NRZ Exp(f)');
title('DSP expérimentale de NRZ en fonction de f');
% pwelch à tester plus tard

%3.1.2.1 Génération du signal x(t) modulé en fréquence
phi0 = rand*2*pi;
phi1 = rand*2*pi;
Cos0 = cos(2*pi*F0*T + phi0);
Cos1 = cos(2*pi*F1*T + phi1);
%Un = eye(size(NRZ));
x_module = (1 - NRZ).*Cos0 + NRZ.*Cos1;

%3.1.2.2 Tracé du signal x(t) en fonction du temps
figure; % figure 3
plot(T,x_module);
ylim([-1.5 1.5]);
xlim([0.05 0.08]);
xlabel('t en s');
ylabel('x(t)');
title('Signal x(t) en fonction de t');

%3.1.2.3 Calcul de la DSP de x(t) en fonction de la DSP de NRZ
% DSP_NRZ_plus_F1 = DSP_rectangulaire(Nb_echantillons,Fe,NRZ,F1);
% DSP_NRZ_plus_F0 = DSP_rectangulaire(Nb_echantillons,Fe,NRZ,F0);
% DSP_NRZ_moins_F1 = DSP_rectangulaire(Nb_echantillons,Fe,NRZ,-F1);
% DSP_NRZ_moins_F0 = DSP_rectangulaire(Nb_echantillons,Fe,NRZ,-F0);
% DSP_X_Exp = (dirac(f-F1) + dirac(f+F1) - dirac(f-F0) - dirac(f-F1))*DSP_NRZ_Exp/4;

%3.1.2.4 Tracé de la DSP de x(t) théorique et expérimentale
DSP_X_Th = DSP_rectangulaire(Nb_echantillons,Fe,x_module,0);

figure; % figure 4
subplot(211);
semilogy(f,DSP_X_Th);
ylim([1e-05 100]);
xlim([-7000 7000]);
xlabel('f en Hz');
ylabel('DSP X Th(f)');
title('DSP X théorique en fonction de f');

subplot(212);
%ylim([1e-04 10]);
%xlim([-6000 6000]);
% semilogy(f,DSP_X_Exp);
% xlabel('f en Hz');
% ylabel('DSP X Exp(f)');
% title('DSP X expérimental en fonction de f');

% S_x = abs(fft(x_module).^2)/length(x_module);
% figure; 
% subplot(411);
% semilogy(S_x);
% subplot(412);
% pwelch(x_module);

%3.2 Canal de transmission à bruit additif, blanc et Gaussien
Px = mean(abs(x_module).^2);
SNR = 10;
sigma = sqrt(Px*10^(-SNR/10));
bruit = sigma*randn(1,Nb_echantillons);

figure; % figure 5 
plot(T,bruit);
ylim([-1.5 1.5]);
xlabel('Temps (secondes)');
ylabel('bruit');
title('Tracé du bruit en fonction du temps');

x_bruite = x_module + transpose(bruit);
figure; % figure 6 
subplot(211);
plot(T,x_module);
ylim([-1.5 1.5]);
xlim([0.05 0.06]);
xlabel('t en s');
ylabel('x(t) sans bruit');
title('Signal x(t) avant ajout du bruit en fonction du temps pour SNR = 10');

subplot(212);
plot(T,x_bruite);
ylim([-1.5 1.5]);
xlim([0.05 0.06]);
xlabel('Temps (secondes)');
ylabel('x bruite');
title('Tracé du signal x bruite en fonction du temps pour SNR = 10');

%3.3 Démodulation par filtrage

%3.3.1 Filtre passe-bas
ordre = 60;
fc = 4000; %Hz
fnorm = fc/Fe;
Tfiltre = (-ordre*Te:Te:ordre*Te);
Ffiltre = (-(2*ordre+1)/2:(2*ordre+1)/2-1)*(Fe/(2*ordre+1));
h_passe_bas = (2*fnorm)*sinc(2*fc*Tfiltre);
y_passe_bas = filter(h_passe_bas,(1),x_module);
H_passe_bas = abs(fftshift(fft(h_passe_bas)));

%3.3.2 Filtre passe-haut
%H_passe_haut = 1 - H_passe_bas;
h_passe_haut = - h_passe_bas;
h_passe_haut(ordre+1) = h_passe_haut(ordre+1) + 1;
y_passe_haut = filter(h_passe_haut,(1),x_module);

%3.3.3 Filtrage ok + manque prise en compte du retard

%3.3.4.1 Tracés réponse impulsionnelle et réponse en fréquence des filtres
figure; % figure 7 
subplot(211);
plot(Tfiltre,h_passe_bas);
title('Réponse impulsionnelle du filtre passe bas');
%ylim([-1.5 1.5]);
%xlim([0.01 0.03]);

subplot(212);
plot(Ffiltre,H_passe_bas);
title('Reponse fréquentielle du filtre passe bas');
%ylim([-1.5 1.5]);
%xlim([0.01 0.03]);

%3.3.4.2 Tracés DSP de x(t) et y(t)
DSP_Y_passe_bas = DSP_rectangulaire(Nb_echantillons,Fe,y_passe_bas,0);
DSP_Y_passe_haut = DSP_rectangulaire(Nb_echantillons,Fe,y_passe_haut,0);

figure; % figure 8 - DSP x(t) et y(t)
subplot(311);
semilogy(f,DSP_X_Th);
ylim([1e-05 100]);
xlim([-7000 7000]);
xlabel('f en Hz');
ylabel('DSP X Th(f)');
title('DSP X théorique en fonction de f');

subplot(312);
semilogy(f,DSP_Y_passe_bas);
ylim([1e-05 100]);
xlim([-7000 7000]);
xlabel('f en Hz');
ylabel('DSP Y passe bas(f)');
title('DSP Y passe bas en fonction de f');

subplot(313);
semilogy(f,DSP_Y_passe_haut);
ylim([1e-05 100]);
xlim([-7000 7000]);
xlabel('f en Hz');
ylabel('DSP Y passe haut(f)');
title('DSP Y passe haut en fonction de f');

% numéro à changer 
figure; % figure 7 - passe bas
subplot(211);
plot(T,x_module);
ylim([-1.5 1.5]);
xlim([0.01 0.03]);

subplot(212);
plot(T,y_passe_bas);
ylim([-1.5 1.5]);
xlim([0.01 0.03]);

figure; % figure 8 - passe haut
subplot(211);
plot(T,x_module);
ylim([-1.5 1.5]);
xlim([0.01 0.03]);

subplot(212);
plot(T,y_passe_haut);
ylim([-1.5 1.5]);
xlim([0.01 0.03]);
% plot(Freq, H_passe_haut);
% plot(intervalle, h_passe_haut);

figure; % figure 9 - passe haut et bas
subplot(311);
plot(T,x_module);
ylim([-1.5 1.5]);
xlim([0.01 0.03]);

subplot(312);
plot(T,y_passe_bas);
ylim([-1.5 1.5]);
xlim([0.01 0.03]);

subplot(313);
plot(T,y_passe_haut);
ylim([-1.5 1.5]);
xlim([0.01 0.03]);
