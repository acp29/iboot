function retval = isparallel

  % Private function file required for ibootci

  % Check we have parallel computing capabilities
  software = ver;
  software = {software.Name};
  if isempty(cell2mat(regexpi(software,'^parallel')))
    retval = false;
  else
    retval = true;
  end

end
