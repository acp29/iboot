% Basic script for local installation
% 

copyfile ('PKG_ADD','PKG_ADD.m');
run (fullfile(pwd,'PKG_ADD.m'));
dirlist = cell(3,1);
dirlist{1} = fullfile (pwd,'inst','helper');
dirlist{2} = fullfile (pwd,'inst','param');
dirlist{3} = fullfile (pwd,'inst','');
n = numel (dirlist);
if isoctave
  % Install for Octave
  octaverc = '~/.octaverc';
  if exist(octaverc,'file')
    [fid, msg] = fopen (octaverc, 'at');
  else
    [fid, msg] = fopen (octaverc, 'wt');
  end
  S = fileread (octaverc);
  comment = '% Load statistics-bootstrap package';
  if isempty(strfind(S,comment))
    fprintf (fid, '\n%s', comment);
  end
  for i=1:n
    if isempty(strfind(S,dirlist{i}))
      fprintf (fid, '\n%s', dirlist{i});
    end
  end
  fclose (fid);
else
  % Assumming install for Matlab instead
  if exist('savepath')
    savepath;
  else
    % backwards compatibility
    path2rc;
  end
end

% Notify user that installation is complete
disp ('The statistics-bootstrap package has been installed at the current location ')

% Clean up
clear dirlist
delete ('PKG_ADD.m');