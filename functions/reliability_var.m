function [v_g v_gs v_gse] = reliability_var(C, )
% Description:
%       variance decomposition on dataset. Say the multi-channel data
%       y_ij follow the following format:
%               y_ij = g + s_i + e_ij
%       i, being the subject number and j, the partition/session number.
%
%       Assuming a) g, s_i, e_ij are mutually independent b) e_ij and s_i
%       are i.i.d, we can estimate the term variances as follows:
%       
%       Across subjects:
%       v_g = E[y_ij, y_kl]
%       Within subject, Across run:
%       v_g + v_s = E[y_ij, y_ik]
%       Within observation:
%       v_g + v_s + v_e = E[y_ij, y_ij]
%
%       To develop estimators for these quantities we replace the 
%       Expectation with the mean over all possible pairings.
%
% INPUT:
%       
