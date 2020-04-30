function y = split_blocks (x, l)

  % Helper function file required for ibootci

  % Calculate data and block dimensions
  n = size(x{1},1);
  nvar = numel(x);

  % Create a matrix of circular, overlapping blocks
  % Ref: Politis and Romano (1991) Technical report No. 370
  y = cell(1,nvar);
  for v = 1:nvar
    y{v} = zeros(n,l);
    temp = cat(1,x{v},x{v}(1:l-1));
    for i = 1:n
      y{v}(i,:) = temp(i:i+l-1);
    end
  end
  y = cell2mat(y);
  y = num2cell(y,1);

end
