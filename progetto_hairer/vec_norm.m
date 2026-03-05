function norm_val = vec_norm(vec, sc)
% VEC_NORM Computes a weighted vector norm.
%
% Syntax:
%   norm_val = vec_norm(vec, sc)
%
% Description:
%   Computes a weighted root-mean-square norm for a vector 'vec',
%   scaled by the components of 'sc'. This is the norm defined in
%   Hairer & Nörsett, formula (4.11), but applied to a single vector
%   instead of a difference of two vectors.
%
% Input:
%   vec - The vector for which to compute the norm.
%   sc  - A vector of scaling factors, same size as 'vec'.
%
% Output:
%   norm_val - The computed scalar norm.

    % Calculate the squared value for each component, weighted by the scale
    val_sq = (vec ./ sc).^2;
    
    % Compute the root-mean-square of the weighted values
    norm_val = sqrt(mean(val_sq));

end
