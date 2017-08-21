function y = fun_sigm(v,varargin)

% this function implements the sigmoid function. Depending on the user
% input, it assumes the values of JuR or user values

    if isempty(varargin)
        'apply default parameters'
        e0=2.5;
        r=560;
        v0=6e-3;
    else
        'apply user parameters'
        e0=varargin{1};
        r=varargin{2};
        v0=varargin{3};
    end
        
    y = 2*e0./(1+exp(r*(v0 - v)));
return;

