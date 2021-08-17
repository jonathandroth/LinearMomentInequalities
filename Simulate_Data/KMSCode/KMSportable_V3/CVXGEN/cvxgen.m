% Creates a local directory called cvxgen_CODENUM, and places a solver called
% csolve in your local directory.

% TODO: make platform independent.

% NOTE: This was updated in 2017-05 to support R2017 of Matlab, using the
% `websave` function. This allows a secure connection to CVXGEN, but will fail
% on older versions of Matlab. Please contact bugs@cvxgen.com if you have
% trouble with this function.

function cvxgen(code)

code = num2str(code);
webopts = weboptions('CertificateFilename', '', 'Timeout', 30);

url = ['https://cvxgen.com/matlab_test/' code];
if ~strcmp(webread(url, webopts), 'success')
  disp(['Failed to retrieve problem ' code '.']);
  return;
else
  disp('Retrieving solver from https://cvxgen.com/...');
end

url = ['https://cvxgen.com/matlab/' code];
dir = 'cvxgen/';

if exist(dir) ~= 7
  status = mkdir(dir);
end

% Download first to cvxgen.zip.
websave('cvxgen.zip', url, webopts);
unzip('cvxgen.zip', '.');

if exist([dir '/solver.c']) == 2
  disp(['Downloaded to ' dir '.']);
else
  disp('Failed to retrieve file. Giving up.');
  return;
end

disp(' ');
disp('cvxgen is in beta. Please note that there is no warranty.');
disp('You may want to check csolve(params) against cvxsolve(params), which uses cvx.');
disp(' ');
disp('Compiling...');

cd(dir);
try
  make_csolve;
  success = 1;
catch exception
  disp('!!! Compilation failed.')
  cd('..');
  success = 0;
  rethrow(exception);
end

if success
  pause(0.05); % Attempt to avoid a certain race condition.
  cd('..');
  pause(0.05);

  pause(0.05); % Attempt to avoid a certain race condition.
  copyfile([dir '/csolve.m*'], '.');
  copyfile([dir '/cvxsolve.m'], '.');

  disp('Success. Type help csolve, or')
  disp(' ');
  disp('  [vars, status] = csolve(params, settings)');
  disp(' ');
end


