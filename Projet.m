clear all
close all
clc
load DonneesBinome1;
% Paramètres
F0 = 6000; % Fréquence bits 0 en Hz
F1 = 2000; % Fréquence bits 1 en Hz
Fe = 48000; % Fréquence d'échantillonage en Hz
Te = 1/Fe; % Période d'échantillonage en secondes
Ts = 1/300; % s/bits, le débit souhaité max est de 300 bits/seconde

Ns = floor(Ts/Te); % échantillons/bits, on prend la partie entière, Ts = Ns*Te
Nb_bits = length(bits); % Nombre de bits du signal émis
%Nb_bits = 100;
Nb_echantillons = Nb_bits*Ns; % Nombre d'échantillons
Null = 2000;

%donnees = randi(0:1, Nb_bits, 1); % Génération de Nb_bits 0 ou 1 de manière aléatoire
donnees = bits;
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
ylim([1e-02 100]);
xlim([-2000 2000]);
xlabel('f en Hz');
ylabel('DSP NRZ Th(f)');
title('DSP théorique de NRZ en fonction de f');

subplot(212);
semilogy(f,DSP_NRZ_Exp);
ylim([1e-6 1e-2]);
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

DSP_X_bruite = DSP_rectangulaire(Nb_echantillons,Fe,x_bruite,0);
figure; % figure 6 bis
semilogy(f,DSP_X_Th);
ylim([1e-01 50]);
xlim([-7000 7000]);
xlabel('f en Hz');
ylabel('DSP X bruite(f)');
title('DSP X bruite en fonction de f');

%3.3 Démodulation par filtrage

%3.3.1 Filtre passe-bas
ordre = 60;
fc = 4000; %Hz
fnorm = fc/Fe;
Tfiltre = (-ordre*Te:Te:ordre*Te);
Ffiltre = (-(2*ordre+1)/2:(2*ordre+1)/2-1)*(Fe/(2*ordre+1));
h_passe_bas = (2*fnorm)*sinc(2*fc*Tfiltre);
y_passe_bas = filter(h_passe_bas,(1),x_bruite);
H_passe_bas = abs(fftshift(fft(h_passe_bas)));

%3.3.2 Filtre passe-haut
H_passe_haut = 1 - H_passe_bas;
h_passe_haut = - h_passe_bas;
h_passe_haut(ordre+1) = h_passe_haut(ordre+1) + 1;
y_passe_haut = filter(h_passe_haut,(1),x_bruite);

%3.3.3 Filtrage ok + manque prise en compte du retard

%3.3.4.1 Tracés réponse impulsionnelle et réponse en fréquence des filtres
figure; % figure 7 
subplot(211);
plot(Tfiltre,h_passe_bas);
title('Réponse impulsionnelle du filtre passe bas');
%ylim([-1.5 1.5]);
xlim([-0.001 0.001]);

subplot(212);
plot(Ffiltre,H_passe_bas);
title('Reponse fréquentielle du filtre passe bas');
%ylim([-1.5 1.5]);
xlim([-9000 9000]);

%3.3.4.2 Tracés DSP de x(t) et réponse fréquentielle du filtre
figure; % figure 8 - DSP x(t) et y(t) passe bas
semilogy(f,DSP_X_Th);
ylim([1e-2 30]);
xlim([-8000 8000]);
xlabel('f en Hz');
ylabel('DSP X Th(f)');
title('DSP X théorique en fonction de f');

hold on;
semilogy(Ffiltre,H_passe_bas);
ylim([1e-02 30]);
xlim([-9000 9000]);

figure; % figure 8 bis - DSP x(t) et y(t) passe haut
semilogy(f,DSP_X_Th);
ylim([1e-2 30]);
xlim([-8000 8000]);
xlabel('f en Hz');
ylabel('DSP X Th(f)');
title('DSP X théorique en fonction de f');

hold on;
semilogy(Ffiltre,H_passe_haut);
ylim([1e-02 30]);
xlim([-9000 9000]);

%3.3.4.3 
figure; % figure 9 - y(t) et DSP_y(t) passe bas
subplot(211);
plot(T,y_passe_bas);
ylim([-1.5 1.5]);
xlim([0.1 0.15]);
xlabel('t en secondes');
ylabel('y(t) signal filtré');
title('Signal y(t) filtré en fonction du temps');

DSP_Y_passe_bas = DSP_rectangulaire(Nb_echantillons,Fe,y_passe_bas,0);
subplot(212);
semilogy(f,DSP_Y_passe_bas);
ylim([1e-2 100]);
xlim([-8000 8000]);
xlabel('f en Hz');
ylabel('DSP Y (f)');
title('DSP Y en fonction de f');

figure; % figure 9 bis - y(t) et DSP_y(t) passe haut
subplot(211);
plot(T,y_passe_haut);
ylim([-1.5 1.5]);
xlim([0.1 0.15]);
xlabel('t en secondes');
ylabel('y(t) signal filtré');
title('Signal y(t) filtré en fonction du temps');

DSP_Y_passe_haut = DSP_rectangulaire(Nb_echantillons,Fe,y_passe_haut,0);
subplot(212);
semilogy(f,DSP_Y_passe_haut);
ylim([1e-2 100]);
xlim([-8000 8000]);
xlabel('f en Hz');
ylabel('DSP Y (f)');
title('DSP Y en fonction de f');


%3.3.5 Détection d'énergie

K = 0.27*Ns;

detect_passe_bas = zeros(Nb_bits,1);
detect_passe_haut = ones(Nb_bits,1);
% On parcourt tous les bits
for i = 1:Nb_bits
    xi_bas = 0;
    xi_haut = 0;
    for j = 1:Ns
        % On calcul la somme des carrés
        xi_bas = xi_bas + y_passe_bas((i-1)*Ns + j)^2;
        xi_haut = xi_haut + y_passe_haut((i-1)*Ns + j)^2;
    end
    % On compare avec le seuil K
    if xi_bas > K
        detect_passe_bas(i) = 1;
    end
    if xi_haut > K
        detect_passe_haut(i) = 0;
    end
end

% Calcul de l'erreur
erreur_passe_bas = sum((transpose(donnees) - detect_passe_bas)).^2;
erreur_passe_haut = sum((transpose(donnees) - detect_passe_haut)).^2;
erreur = mean(erreur_passe_bas, erreur_passe_haut)

% Reconstitution image

suite_binaire_reconstruite = detect_passe_bas;

pcode reconstitution_image;
reconstitution_image(suite_binaire_reconstruite) ;
which reconstitution_image;