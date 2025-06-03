function gabor = generateGabor(rfRad, theta)
    

    sigma = sqrt(rfRad); 
    theta_rad = deg2rad(theta); 
    
    
    s = 2 * rfRad + 1; % Scale of filter
    [X, Y] = meshgrid(-rfRad:rfRad, -rfRad:rfRad);
    
    % Rotation
    x_theta =  X * cos(theta_rad) + Y * sin(theta_rad);
    y_theta = -X * sin(theta_rad) + Y * cos(theta_rad);
    
    
    gamma = 0.5; % Eccentricity
    psi = 0; % Phase

    gaussian = exp(-(x_theta.^2 + (gamma^2) * y_theta.^2) / (2 * sigma^2));
    cosine = cos(2 * pi * x_theta / rfRad + psi);
    
    % Gabor Filter
    gabor = gaussian .* cosine;
    
    %Normalization
    gabor = gabor / sum(abs(gabor()),'all'); 
end
