function [dXi_de, dAlpha_da] = dXiAlpha_deda(e, a, omega)
    % Input:
    % e: Nr x 3 matrix
    % a: Nr x 4 matrix
    % omega: scalar or vector

    % Number of rows
    Nr = size(e, 1);

    % Compute partial derivatives for xi_omega
    % xi_omega = e(:,1) ./ (1 + e(:,2) .* omega^2).^e(:,3);
    xi_omega_denom = (1 + e(:,2) * (omega^2)); % Adjusted for omega^2
    xi_omega_denom_pow = xi_omega_denom.^e(:,3);
    
    dXi_de = zeros(Nr, 3);
    dXi_de(:,1) = 1 ./ xi_omega_denom_pow;
    dXi_de(:,2) = -e(:,1) .* e(:,3) .* (omega^2) ./ (xi_omega_denom_pow .* xi_omega_denom); % Adjusted for omega^2
    dXi_de(:,3) = -e(:,1) .* log(xi_omega_denom) ./ xi_omega_denom_pow;

    % Compute partial derivatives for alpha_omega
    % alpha_omega = a(:,1) ./ (1 + a(:,2) .* ((omega - a(:,4))^2)).^a(:,3);
    alpha_omega_denom = (1 + a(:,2) .* ((omega - a(:,4)).^2)); % Adjusted for (omega - a(:,4))^2
    alpha_omega_denom_pow = alpha_omega_denom.^a(:,3);
    
    dAlpha_da = zeros(Nr, 4);
    dAlpha_da(:,1) = 1 ./ alpha_omega_denom_pow;
    dAlpha_da(:,2) = -a(:,1) .* a(:,3) .* ((omega - a(:,4)).^2) ./ (alpha_omega_denom_pow .* alpha_omega_denom); % Adjusted for (omega - a(:,4))^2
    dAlpha_da(:,3) = -a(:,1) .* log(alpha_omega_denom) ./ alpha_omega_denom_pow;
    dAlpha_da(:,4) = 2 * a(:,1) .* a(:,3) .* a(:,2) .* (omega - a(:,4)) ./ (alpha_omega_denom_pow .* alpha_omega_denom); % Derivative adjusted with factor 2(omega-a(:,4))
end
