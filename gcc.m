
function [zc]=gcc(z1, z2, flag)
% Code taken from https://github.com/javiribera/TDE-and-whale-localization/tree/master/algorithms_TDE
% |gcc(z1, z2, flag)| computes the generalized cross correlation (GCC)
% between signals z1 and z2, from FFT/IFFT, as specified in (Knapp & Carter
% 1976).  
% |[flag]| makes it possible to choose the type of cross-correlation:
%  Standard cross correlation if |flag='cc'|;
%  Phase transform:if |flag='phat'|;
%  Smoothed coherence transform:if |flag='scot'|
%
% C. H. Knapp and G.C. Carter, "The Generalized Correlation Method for
% Estimation of Time Delay", IEEE, Trans. on ASSP, No4(24), Aug. 1976

% Original Source from Thierry Dutoit, Ferran Marquès, 
%   "Applied Signal Processing-A MATLAB-Based Proof of Concept"
%   Springer:New-York, 2009

% make sure we work with rows
z1 =z1(:)';
z2 =z2(:)';

M = min(length(z1),length(z2)); % length of signal
NFFT = 2*M-1; 
z1 = z1(1:M);
z1 = z1 - mean(z1);
z2 = z2(1:M); 
z2 = z2 - mean(z2);
Z1 = fft(z1,NFFT);
Z2 = fft(z2,NFFT);

% cross correlation
Phi_z1z2 = Z1.*conj(Z2);
if(strcmp(flag,'cc'))
   phi_z1z2 = ifft(Phi_z1z2);
   zc = [phi_z1z2(NFFT-M+2:NFFT) phi_z1z2(1:M)];  % re-arrange the vector
elseif(strcmp(flag,'phat'))
   phi_z1z2 = ifft(Phi_z1z2 ./ max(abs(Phi_z1z2),eps));
   zc = [phi_z1z2(NFFT-M+2:NFFT) phi_z1z2(1:M)];  % re-arrange the vector

elseif(strcmp(flag,'scot'))
    % Smoothed Coherence Transform (SCOT)  
    % this processor exhibits the same spreading as the Roth processor. 

    W = 1./( (Z1.*Z2).^0.5 );
    phi_z1z2 = ifft(Phi_z1z2 .* W);
    zc = [phi_z1z2(NFFT-M+2:NFFT) phi_z1z2(1:M)];  % re-arrange the vector  

else
   disp('Invalid value for flag');
end

zc = real(zc); %just in case (this way we can visualize gcc scot)

