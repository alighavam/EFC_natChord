function [n,d] = nSphere_count_samples(nat_dist,chord,r,varargin)
% Description:
%       Creates an K dimensional (K being number of EMG channels) n-sphere 
%       with radius r centering chord. Then counts the number of samples
%       from nat_dist that falls in the n-sphere. (This description is the
%       intuition of what the function does. But in fact, an n-sphere is
%       not actually created and the samples are counted by using a
%       distance measure.
%
%       distance = (x_i - u)' sigma^-1 (x_i - u)
%       x_i: a point in nat_dist. 
%       u: chord pattern
%       sigma: depends on the 'd_type'. Refer to INPUT section.
%       
% INPUT:
%       nat_dist: natural EMG distribution. Rows are observations and
%       Columns are EMG channels. 
%       
%       chord: chord EMG pattern. It's a column vector of chord
%       EMGs. Must be an K dimensional vector, K being number of channels.
%
%       r: radius of the n-sphere. 
%
%   varargin:
%       'd_type' -> 'Euclidean', 'project_to_nSphere', 'oval' 
%
%       type of the distance to use to count the number of samples.
%       Intuition of different distance types is the shape of the n-sphere.
%
%       Euclidean: sigma = I (N by N) -> an n-sphere with equal radius in
%       all dimensions.
%
%       oval: sigma = I + lambda * (u*u') -> an n-sphere with the shape of
%
%       oval in some dimensions. The shape depends on the covariance matrix
%       of the chord pattern.
%
%       project_to_nSphere: project all points to the n-sphere and calculate
%       Euclidean distance on the n-sphere.
%       
%       'lambda' -> parameter for 'oval' d_type we need to assign a lambda.

% handling input arguments:
d_type = 'Euclidean';
lambda = [];
vararginoptions(varargin,{'d_type','lambda'})

if (~isvector(chord))
    error('nSphere_count_samples: chord must be a vector')
end
if (~iscolumn(chord))
    chord = chord';
end

switch d_type
    case 'Euclidean'
        lambda = 0;

    case 'oval'
        if (isempty(lambda))
            warning(['nSphere_count_samples: When using oval distance ' ...
                     'option, you must input a lambda. Setting lambda to 1' ...
                     ' and calculating the oval distance.'])
            lambda = 1;
        end

    case 'project_to_nSphere'
        lambda = 20000;

    otherwise
        error('nSphere_count_samples: Distance %s does not exist.',d_type)
end

% distance container:
d = zeros(size(nat_dist,1),1);

% variance covariance matrix of the chord pattern:
cov_chord = chord * chord';

% distance weights:
sigma = eye(size(cov_chord)) + lambda * (cov_chord);

% initialize number of samples:
n = 0;

% looping through points and calculating distances:
for i = 1:size(nat_dist,1)
    % sample in natural distribution:
    x = nat_dist(i,:)';

    % squared distance:
    d(i) = (x - chord)' * sigma^-1 * (x - chord);

    % disp(d(i))

    % counting if sample lies in the n-sphere:
    n = n + double(d(i) < r^2); 
end







