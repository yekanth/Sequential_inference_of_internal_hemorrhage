function HCT_with_noise = add_noise(HCT,noise,noise_type)
    
    % Function to add various types of noises to the simulated measurements
    if noise_type
        HCT_with_noise = HCT + noise .* randn(size(HCT));
    else 
        HCT_with_noise = HCT + noise .* (-1 + rand(size(HCT)));
    end

end