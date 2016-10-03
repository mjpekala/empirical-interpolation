function setup()
% SETUP  Update Matlab search path.
%
% References:
%   [GPML] http://www.gaussianprocess.org/gpml/code/matlab/doc/

% mjp, oct 2016

% update appsDir as needed for the local system
here = fileparts(mfilename('fullpath'));
app_dir = fullfile('/Users', 'pekalmj1', 'Apps');


if exist('covADD') ~= 2
    gpml_dir = fullfile(app_dir, 'gpml-matlab-v3.6-2015-07-07');
    run(fullfile(gpml_dir, 'startup.m'));
    fprintf('[%s]: using "%s" for GPML\n', mfilename, gpml_dir);
end
