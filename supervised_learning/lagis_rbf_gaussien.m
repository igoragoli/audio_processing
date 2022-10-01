%%%%%%%%%%%%%%%%%%%%
% SVM - noyau RBF gaussien
% Code Lagis et CentraleSupélec
%
% Arguments :
% -----------
% 
% - "donnee1" et "donnee2" sont des matrices. Deux cas se présentent :
%    1) Elles ont même dimension :
%    * "N" lignes, où "N" est le nombre de paramètres (features)
%      extraits, c'est-à-dire encore "ndimensions"
%    * "m" colonnes, où "m" est le nombre d'observations extraites
%    2) L'une est un vecteur d'une colonne et l'autre est comme
%       indiqué en 1)
% - "sig2" est calculé séparément ("lagis_sig2.m" pour RBF gaussien)
%    mais ce n'est pas encore fourni en tant que fonction standard,
%    alors "sig2" doit encore être entré en paramètre
%
% Stéphane Rossignol -  14/02/2006 ; 2021 pour CentraleSupélec
%
%%%%%%%%%%%%%%%%%%%%

function [valeurs] = lagis_rbf_gaussien (donnee1, donnee2, sig2)

Ndim = size(donnee1,1);

normes = 0;
for ii=1:Ndim
   normes = normes + (donnee1(ii,:)-donnee2(ii,:)).^2;
end;

valeurs = exp ( -normes/sig2 );
