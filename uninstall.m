% Basic uninstall script 
% 


dirlist = fullfile (pwd,'inst',{"", "helper", "param"});
N = cellfun('numel',dirlist);
[S, I] = sort (N, 'descend');
sortedlist = dirlist(I);
n = numel(dirlist);
if isoctave
  octaverc = '~/.octaverc';
  if exist(octaverc,'file')
    S = fileread(octaverc);
    [fid, msg] = fopen (octaverc, 'w');
  else
    error('~/.octaverc does not exist');
  end
  comment = '\n% Load statistics-bootstrap package';
  S = regexprep(S,comment,'');
  for i=1:n
    S = regexprep(S,sprintf('\n%s',sortedlist{i}),'');
  end
end
fwrite (fid,S,'char');
fclose (fid);

disp(sprintf('statistics-bootstrap package has been uninstalled from: %s ',dirlist{1}))

run('PKG_DEL');