function L(k, q, λ)
    a, b, n = (k - 0.5) * q, (k + 0.5) * q, 20
    delta = (b - a) / n
    xs = a .+ (0:n-1) * delta
    return (λ/2) * sum(exp.(-λ .* abs.(xs))) * delta
    return
end

function likelihood_given_q(t, q, σ_2, ζ, λ)
    # Need to calculate the k range
    # This means the min value of k is such that t - kq = ζ => k = round_towards_zero ((t - ζ)/q)
    # And the max value of k is such that t - kq = -ζ => k = round_towards_zero((t + ζ)/q)
    bound_1_k = (t - ζ) / q
    bound_1_k = bound_1_k < 0 ? ceil(bound_1_k) : floor(bound_1_k)

    bound_2_k = (t + ζ) / q
    bound_2_k = bound_2_k < 0 ? ceil(bound_2_k) : floor(bound_2_k)

    k = Int64(min(bound_1_k, bound_2_k)):Int64(max(bound_1_k, bound_2_k))

    if length(k) == 0
        return 0
    else
        return sum(exp.((-(t .- q .* k) .^ 2) ./ (2 * σ_2)) .* L.(k, q, λ))
    end
end

function refined_grayscale_estimation(m, n, coefficients, potential_q_values)
    Di = coefficients[:, m, n]
    λ = length(Di) / sum(abs.(Di))
    σ_2 = 0.8
    ζ = 6
    best_q = 0
    highest_likelihood = -Inf
    for q in potential_q_values[:, m, n]
        curLikelihood = sum(log.(likelihood_given_q.(Di, q, σ_2, ζ, λ)))
        if curLikelihood > highest_likelihood
            highest_likelihood = curLikelihood
            best_q = q
        end
    end
    return best_q
end

function calculate_potential_q_values(quantisation_table)
    potential_q_values = Int64.(zeros(100,8,8))
    for qual in 1:100
        potential_q_values[qual, :, :] = Int64.(round.(JPEG_quality_multiplier.(qual) .* quantisation_table))
    end
    replace!(potential_q_values, 0 => 1)
    return potential_q_values
end

function calculate_refined_mle_matrix(coefficients, potential_q_values)
    refined_mle_matrix = zeros(8,8)
    for m in 1:8
        for n in 1:8
            refined_mle_matrix[m,n] = refined_grayscale_estimation(m, n, coefficients[:, :, :], potential_q_values)
        end
    end
    return refined_mle_matrix
end

function remove_smoothing(smoothed_image)
    # This function attempts to remove the smoothing which is done by DJEG
    # The smoothing is done by convolving with 
    # 1 2 1
    # 2 4 2 all times 1/16
    # 1 2 1

    # First buildup a zero padded version of the image
    zero_padded_image = zeros(size(smoothed_image) .+ 4)
    zero_padded_image[3:end-2, 3:end-2] = smoothed_image
    zero_padded_image[3:end-2,2] = 3/4 .* zero_padded_image[3:end-2,3]
    zero_padded_image[3:end-2,1] = 1/4 .* zero_padded_image[3:end-2,3]
    zero_padded_image[3:end-2,end-1] = 3/4 .* zero_padded_image[3:end-2,end-2]
    zero_padded_image[3:end-2,end] = 1/4 .* zero_padded_image[3:end-2,end-2]
    zero_padded_image[2,:] = 3/4 .* zero_padded_image[3,:]
    zero_padded_image[1,:] = 1/4 .* zero_padded_image[3,:]
    zero_padded_image[end-1,:] = 3/4 .* zero_padded_image[end-2,:]
    zero_padded_image[end,:] = 1/4 .* zero_padded_image[end-2,:]

    # Make the smoothing filter, the impulse response
    h = zeros(size(zero_padded_image))
    h[1, 1] = 4
    h[1, 2] = 2
    h[2, 1] = 2
    h[2, 2] = 1
    h[end, end] = 1
    h[1, end] = 2
    h[end, 1] = 2
    h[end, 2] = 1
    h[2, end] = 1
    h ./= 16

    # Now attempt to undo the convolution
    Y = fft(zero_padded_image)
    H = fft(h)

    len1 = size(zero_padded_image)[1]
    len2 = size(zero_padded_image)[2]

    PSD = abs.(Y) .^ 2
    sigma = 0.33
    wiener = Y .* ((conj.(H) .* PSD) ./ ((abs.(H .^ 2) .* PSD) .+ (len1 * len2) * sigma^2))

    PSD = abs.(wiener .^ 2) .+ 0.1 * sigma^2
    wiener = Y .* ((conj.(H) .* PSD) ./ ((abs.(H .^ 2) .* PSD) .+ (len1 * len2) * sigma^2))

    estimate = real(ifft(wiener))
    recovered = estimate[3:end-2, 3:end-2]
    recovered1 = conv(recovered, ones(2,2) ./ 4)

    recovered2 = recovered1[2:2:end, 2:2:end]
    recovered3 = zeros(size(estimate[3:end-2, 3:end-2]))
    recovered3[1:2:end,1:2:end] = recovered2;
    recovered3[1:2:end,2:2:end] = recovered2;
    recovered3[2:2:end,1:2:end] = recovered2;
    recovered3[2:2:end,2:2:end] = recovered2;
    origZPimage = zeros(size(recovered3) .+ 4);
    origZPimage[3:end-2,3:end-2] = recovered3;
    origZPimage[3:end-2,2] = origZPimage[3:end-2,3];
    origZPimage[3:end-2,end-1] = origZPimage[3:end-2,end-2];
    origZPimage[2,:] = origZPimage[3,:];
    origZPimage[end-1,:] = origZPimage[end-2,:];

    PSD = abs.(fft(origZPimage) .^ 2) .+ 0.1 * sigma^2
    wiener = Y .* ((conj.(H) .* PSD) ./ ((abs.(H .^ 2) .* PSD) .+ (len1 * len2) * sigma^2))
    estimate = real(ifft(wiener))
    recovered = estimate[3:end-2, 3:end-2]
    recovered1 = conv(recovered, ones(2,2) ./ 4)
    return recovered1[2:2:end, 2:2:end]
end