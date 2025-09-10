function [dXi_de, dAlpha_da] = dXiAlpha_deda(e, a, omega)
% Input:
    % e: Nr x 3 matrix
    % a: Nr x 4 matrix
    % omega: scalar or vector
    
    % Number of rows
    Nr = size(e, 1);
    %% --- Xi derivatives ---
    xi_omega_denom = (1 + e(:,2) * (omega^2));
    xi_omega_denom_pow = xi_omega_denom.^e(:,3);
    
    dXi_de = zeros(Nr, 3);
    dXi_de(:,1) = 1 ./ xi_omega_denom_pow;
    dXi_de(:,2) = -e(:,1) .* e(:,3) .* (omega^2) ./ (xi_omega_denom_pow .* xi_omega_denom);
    dXi_de(:,3) = -e(:,1) .* log(xi_omega_denom) ./ xi_omega_denom_pow;
    %% --- Alpha derivatives ---
    % alpha_omega = 1/2*( term_pos + term_neg )
    denom_pos = (1 + a(:,2) .* ( (omega - a(:,4)).^2 ));
    denom_neg = (1 + a(:,2) .* ( (-omega - a(:,4)).^2 ));
    denom_pos_pow = denom_pos.^a(:,3);
    denom_neg_pow = denom_neg.^a(:,3);
    dAlpha_da = zeros(Nr, 4);
    % derivative wrt a1
    dAlpha_da(:,1) = 0.5 * ( 1 ./ denom_pos_pow + 1 ./ denom_neg_pow );
    % derivative wrt a2
    dAlpha_da(:,2) = 0.5 * ( ...
        -a(:,1) .* a(:,3) .* ((omega - a(:,4)).^2) ./ (denom_pos_pow .* denom_pos) ...
        -a(:,1) .* a(:,3) .* ((-omega - a(:,4)).^2) ./ (denom_neg_pow .* denom_neg) );
    % derivative wrt a3
    dAlpha_da(:,3) = 0.5 * ( ...
        -a(:,1) .* log(denom_pos) ./ denom_pos_pow ...
        -a(:,1) .* log(denom_neg) ./ denom_neg_pow );
    % derivative wrt a4
    dAlpha_da(:,4) = 0.5 * ( ...
        2 * a(:,1) .* a(:,3) .* a(:,2) .* (omega - a(:,4)) ./ (denom_pos_pow .* denom_pos) ...
      + 2 * a(:,1) .* a(:,3) .* a(:,2) .* (-omega - a(:,4)) ./ (denom_neg_pow .* denom_neg) );
end