function d = get_d_from_natural(pattern,nat_dist,varargin)
% Description:
%       Calculates distance of pattern from every chord in nat_dist
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
%       pattern: EMG pattern. It's a column vector of EMGs. Must be a 
%       K by 1 vector, K being number of EMG channels.
%
%   varargin:
%       'd_type' -> 'Euclidean', 'project_to_nSphere', 'oval' 
%
%       type of the distance to use to count the number of samples.
%       Intuition of different distance types is the projection of the
%       natural distribution samples to the n-spher.
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

if (~isvector(pattern))
    error('get_d_from_natural: pattern must be a vector')
end
if (~iscolumn(pattern))
    pattern = pattern';
end

switch d_type
    case 'Euclidean'
        lambda = 0;

    case 'oval'
        if (isempty(lambda))
            warning(['get_d_from_natural: When using oval distance ' ...
                     'option, you must input a lambda. Setting lambda to 1' ...
                     ' and calculating the oval distance.'])
            lambda = 1;
        end

    case 'project_to_nSphere'
        lambda = 20000;

    otherwise
        error('get_d_from_natural: Distance %s does not exist.',d_type)
end

% distance container:
d = zeros(size(nat_dist,1),1);

% variance covariance matrix of the chord pattern:
cov_pattern = pattern * pattern';

% distance weights:
sigma = eye(size(cov_pattern)) + lambda * (cov_pattern);

% looping through points and calculating distances:
for i = 1:size(nat_dist,1)
    % sample in natural distribution:
    x = nat_dist(i,:)';

    % squared distance:
    d(i) = (x - pattern)' * sigma^-1 * (x - pattern);
end

d = sqrt(d);
d = sort(d);




