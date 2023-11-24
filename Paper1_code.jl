function get_DCT_coefficients_excluding_truncated_uniform(image_folder, image_name)
    base_image = Gray.(FileIO.load(image_folder * image_name * ".png"))
    base_image_integers = Int64.(round.(Float64.(base_image) .* 255))
    dct_coefficients = zeros( 1 + (length(base_image) ÷ 64), 8, 8)
    block_index = 1
    for i in 1:8:axes(base_image_integers)[1][end]
        for j in 1:8:axes(base_image_integers)[2][end]
            block = Float64.(base_image_integers[i:i+7, j:j+7])
            if maximum(block) == 255 || minimum(block) == 0 || maximum(block) == minimum(block)
            else
                block_coefficients = FFTW.dct(block .- 128)
                dct_coefficients[block_index, :, :] = round.(block_coefficients)
                block_index += 1
            end
        end
    end
    return dct_coefficients[1:block_index, :, :]
end

function D(m)
    if m == 0 || m == 4
        return 2
    elseif m == 2 || m == 6
        return 2 * cos(pi / 4)
    else
        return 2 * cos(pi / 4) * cos(pi / 8)
    end
end

B(m,n) = D(m) * D(n)

function G(x, B)
    if abs(x) > B
        return 0
    else
        return exp(-6 * (-x)^2) # Ignoring the normalisation constant because 
    end
end

function integral_of_sum_G(Ys, r, q, b)
    # For the time being assume that k is in the range 0..255

    max_k = Int64(ceil(maximum(Ys + b) / q))
    min_k = Int64(floor(min(Ys - b) / q))
    k = -min_k:max_k
    sum_G(x) = sum(G.(x .- r*q .- k.*q, b))

    # Approximate the integral by using a left reimann sum
    a, b, n = Ys - 0.5, Ys + 0.5, 20
    delta = (b - a) / n
    xs = a .+ (0:n-1) * delta
    return sum(sum_G.(xs)) * delta
end

log_likelihood(Y_prime, r, q, b, N) = sum(log.(integral_of_sum_G.(Y_prime, r, q, b) .+ 0.000000001)) + N * log(q)

function mle_q(m, n, all_coefficients) 
    b = B(m, n)
    Y_prime = all_coefficients[:, m, n]
    N = length(Y_prime)
    if length(abs.(Y_prime[abs.(Y_prime) .> b])) == 0 
        return -1
    end
    Q = Int64.(mode(abs.(Y_prime[abs.(Y_prime) .> b])))
    q = [1,2,3,5,Q-1,Q,Q+1]
    max_q_value = 0
    max_val = -Inf
    # q are the values we are trying for quantisation
    # s is iterating through all of the blocks
    # k is like the number of peaks out from the centre we are because each peak is at k times the quantisation value
    # r is round(Y* / q)
    # N is the number of blocks we are approximating over
    for cur_q in q
        # println(cur_q)
        r = round.(Y_prime / cur_q)
        val = log_likelihood(Y_prime, r, cur_q, b, N)
        if val > max_val
            max_val = val
            max_q_value = cur_q
        end
    end
    return max_q_value
end

function likelihood_of_quality(Q, all_coefficients)
    total_likelihood = 0
    quantitsation_table = JPEG_quality_multiplier(Q) .* standard_luminance_table_q50
    for m in 1:8
        for n in 1:8
            if m == 1 && n == 1
            else
                # This is making sure the DC component is excluded from the calculation
                b = B(m, n)
                Y_prime = all_coefficients[:, m, n]
                N = length(Y_prime)
                cur_q = quantitsation_table[m, n]
                r = round.(Y_prime / cur_q)
                total_likelihood += log_likelihood(Y_prime, r, cur_q, b, N)
            end
        end
    end
    return total_likelihood
end

function mle_quality(all_coefficients)
    max_q = 0
    max_likelihood = -Inf
    for q in 1:99
        likelihood = likelihood_of_quality(q, all_coefficients)
        if likelihood > max_likelihood
            max_likelihood = likelihood
            max_q = q
        end
    end
    return max_q
end

function get_mle_quant_table(coefficients)
    mle_quant_table = zeros(8,8)
    for m in 1:8
        for n in 1:8
            mle_quant_table[m, n] = mle_q(m, n, coefficients[:, :, :])
        end
    end
    return mle_quant_table
end
