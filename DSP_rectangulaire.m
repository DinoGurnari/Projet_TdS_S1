
function [DSP_signal,f] = DSP_rectangulaire(Nb_points,Fe,signal,Translation)

DSP_signal = periodogram(signal,[],Nb_points,'centered');
f = (-Nb_points/2:Nb_points/2-1)*(Fe/Nb_points) + Translation;
end