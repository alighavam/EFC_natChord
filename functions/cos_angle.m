function gamma = cos_angle(A,B)
% Description:
%       A and B are vertical vectors or matrices of size NxT. 
%
%       If A and B are vectors, gamma equals to the cosine angle between A
%       and B.
%
%       If A and B are matrices, gamma equals to the cosine angle between
%       the corresponding columns of A and B. 
%       Like col 1 of A and col 1 of B. 
%       
% Output:
%       gamma is a 1xT vector where T depends on the size of A and B.


gamma = dot(A, B)./(vecnorm(A).*vecnorm(B));