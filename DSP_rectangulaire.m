
function [DSP_signal,f] = DSP_rectangulaire(Nb_points,Fe,signal)

DSP_signal = periodogram(signal,[],Nb_points,'centered');
%N = Nb_points;
f = (-Nb_points/2:Nb_points/2-1)*(Fe/Nb_points);
end