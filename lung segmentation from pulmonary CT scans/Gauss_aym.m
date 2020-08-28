
function [fy]=Gauss_aym(gray,mean,segma)

fy=((1./sqrt(2*pi.*segma)))*exp((-(gray-mean)^2)/(2*segma));
end
