% Basic script for local installation
% 

run ('PKG_ADD');


dirlist = fullfile (pwd,'inst',{"", "helper", "param"});
n = numel (dirlist);
if isoctave
  octaverc = '~/.octaverc';
  if exist(octaverc,'file')
    [fid, msg] = fopen (octaverc, 'a');
  else
    [fid, msg] = fopen (octaverc, 'w');
  end
  S = fileread (octaverc);
  comment = '% Load statistics-bootstrap package';
  if isempty(regexp(S,comment,'once'))
    fprintf (fid, '\n%s', comment)
  end
  for i=1:n
    if isempty(regexp(S,dirlist{i},'once'))
      fprintf (fid, '\n%s', dirlist{i})
    end
  end
end
fclose (fid);

disp (sprintf ('statistics-bootstrap package has been installed locally at: %s ',dirlist{1}))