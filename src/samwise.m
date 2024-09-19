function [geo, mat, bcd, sol, opt, res, msh, glb] = samwise(varargin)
%% SAMWISE - Semi-Analytical Modeling of Waves In Structural Elements
% This function checks the inputs and calls the main routine, which
% includes the typical finite element computations.
% The input can be in arbitrary order. One of the inputs must be an object
% of a subclass of "geometry". Everything else will be set to default values if
% undefined. Solutions will always be computed if a specific solver is
% provided.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% check input
% check which objects are provided, and set default values to the others
[geo, mat, bcd, sol, opt, res] = checkInput(varargin);

%% call main routine
[geo, mat, bcd, sol, opt, res, msh, glb] = samwiseMain(geo, mat, bcd, sol, opt, res);

end
