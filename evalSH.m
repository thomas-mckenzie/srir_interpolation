function [Ysh] = evalSH(maxOrder, dirs_sph)
% Evaluate real-valued Spherical Harmonics 
% maxOrder ... maximal Ambisonics order
% dirs_sph ... directions in sperical coordiates: 
%              [azi, ele], angles in radiant 
%              dim(dirs_sph) = Q x 2
% output:  Ysh ... matrix of spherical harmonics 
%          dim(Ysh) = Q x (N + 1)^2

numDirs = size(dirs_sph, 1);
assert(size(dirs_sph, 2)== 2, ...
    'evalSH needs D x 2 matrix with columns [azimuth, elevation]');

phi = dirs_sph(:, 1); % azimuth
theta = dirs_sph(:, 2); % elevation

% initialize SH matrix
Ysh = zeros(numDirs, (maxOrder+1)^2);

% compute SHs for every order and sort them into the matrix
for n = 0:maxOrder
    
    % Azimuthal part
    Yazi = [ sqrt(2) * sin(phi * (n:-1:1)), ...
        ones(numDirs, 1), sqrt(2) *  cos( phi * (1:n))];
    
    % Zenithal part
    Yzen = legendre(n, cos(pi / 2 - theta))';
    csPhase = ones(numDirs,1) * (-1).^(-n:n) ;
    
    % normalization
    normlz = sqrt( (2*n+1)*factorial(n-abs(-n:n)) ./ ...
        (4*pi*factorial(n+abs(-n:n))))';
    
    % put together
    idxSh =  n^2 + n + (-n:n) + 1;
    Ysh(:, idxSh) = (Yazi .* csPhase .* [Yzen(:, end:-1:1), Yzen(:, 2:end)] ) .* normlz';
end

end



