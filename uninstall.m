% Basic uninstall script 
% 

copyfile ('PKG_ADD','PKG_ADD.m'); % copy to m-file for backwards compatibility
copyfile ('PKG_DEL','PKG_DEL.m'); % copy to m-file for backwards compatibility
run (fullfile(pwd,'PKG_ADD.m'));

% Check if running in Octave (else assume Matlab)
info = ver; 
isoctave = any (ismember ({info.Name}, 'Octave'));

if isoctave
  % Uninstall for Octave
  dirlist = cell(2,1); % dir list needs to be in decreasing order of length
  dirlist{1} = fullfile (pwd,'inst','helper');
  dirlist{2} = fullfile (pwd,'inst','');
  n = numel (dirlist);
  octaverc = '~/.octaverc';
  if exist(octaverc,'file')
    [fid, msg] = fopen (octaverc, 'r+t');
    S = (fread (fid, '*char')).';
    fclose(fid);
    [fid, msg] = fopen (octaverc, 'wt');
  else
    error('~/.octaverc does not exist');
  end
  comment = regexptranslate ('escape', '% Load iboot package');
  S = regexprep(S,['\r\n\r\n',comment],'');
  for i=1:n
    S = regexprep(S,strcat('\r\n',...
                  regexptranslate ('escape', strcat('addpath (''',dirlist{i},''', ''-end'')'))),'');
  end
  fseek (fid, 0);
  fputs (fid, S);
  fclose (fid);
  % Unload package from current session
  run (fullfile(pwd,'PKG_DEL.m'));
else
  % Assumming uninstall for Matlab instead
  run (fullfile(pwd,'PKG_DEL'));
  if exist('savepath')
    savepath
  else
    % backwards compatibility
    path2rc;
  end
end


% Notify user that uninstall is complete
disp ('This iboot package has been uninstalled from this location')

% Clean up
clear info isoctave dirlist S comment i ii octaverc fid n msg
delete ('PKG_ADD.m');
delete ('PKG_DEL.m');
