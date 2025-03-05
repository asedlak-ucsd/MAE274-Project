function psys = permute_ss(sys, idx)
    % Permute the order of equations in a state-space model so that all
    % states in idx are the final states
    n = length(sys.A);
    states = 1:n;
    states = [setdiff(states, idx) idx];
    % Generate a permutation matrix P such that P * A will permute the rows 
    % according to the idx array
    P = zeros(n);
    
    for i = 1:n
      P(i, states(i)) = 1;
    end
    
    psys = ss2ss(sys, P);
  end