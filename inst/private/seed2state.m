function state = seed2state(seed)

  % Convert seed to state of the Mersenne Twister random number generator
  state = uint32(zeros(625,1));
  state(1) = seed;
  for i = 1:623
    state(i+1) = uint32(bitand(...
                   (uint64(1812433253) * uint64(bitxor(state(i),bitshift(state(i),-30))) + i),...
                    uint64(intmax('uint32'))));
  end
  if isoctave
    state(end) = 1;
  else
    state(end) = 624;
  end

end
