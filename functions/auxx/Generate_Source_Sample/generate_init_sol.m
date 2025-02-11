function  [x] = generate_init_sol(x_mean,x_std,R)

x_new = x_mean + x_std.*randn(size(x_mean));
x_new=max(x_new,0);
[e,a,s2] = x2v(x_new);

% Identify non-zero entries in a and e
    nonzero_positions_a = a(:, 1) ~= 0;
    nonzero_positions_e = e(:, 1) ~= 0;

%     % Apply conditional modifications based on non-zero positions in a
    if any(nonzero_positions_a)
         a(nonzero_positions_a, 4) = max(min(a(nonzero_positions_a, 4), 13), 7);
         a(nonzero_positions_a, 2) = max(min(a(nonzero_positions_a, 2), 0.009), 0.0009);
         a(nonzero_positions_a, 3) = max(min(a(nonzero_positions_a, 3), 40), 20);
        % a(nonzero_positions_a, 1) = min(a(nonzero_positions_a,1),e(nonzero_positions_a,1)/5);
     end
% 
%     % Apply conditional modifications based on non-zero positions in e
    if any(nonzero_positions_e)
        e(nonzero_positions_e, 2) = max(min(e(nonzero_positions_e, 2), 0.009), 0.0009);
        e(nonzero_positions_e, 3) = max(min(e(nonzero_positions_e, 3), 3), 1.5);
    end

    % Identify zero entries in a and e
    zero_positions_a = a(:, 1) == 0;
    zero_positions_e = e(:, 1) == 0;

    % Set columns 2 to 4 to zero for zero positions in a
    if any(zero_positions_a)
        a(zero_positions_a, 2:4) = 0;
    end

    % Set columns 2 and 3 to zero for zero positions in e
    if any(zero_positions_e)
        e(zero_positions_e, 2:3) = 0;
    end

    % Identify zero entries in a and e
    zero_positions_a2 = a(:, 4) == 0;

    % Set columns 2 to 4 to zero for zero positions in a
    if any(zero_positions_a2)
        a(zero_positions_a2, 1) = 0;
    end
    
   e = R*e;
   a = R*a;
   s2 =s2;
    % Reassemble the parts e, a, s2 back into p using v2x function
    x = v2x(e, a, s2);
end

