
jpeg_coefficient_numbers = [1  2  6  7  15 16 28 29;
 3  5  8  14 17 27 30 43;
 4  9  13 18 26 31 42 44;
 10 12 19 25 32 41 45 54;
 11 20 24 33 40 46 53 55;
 21 23 34 39 47 52 56 61;
 23 35 38 48 51 57 60 62;
 36 37 49 50 58 59 63 64;]

standard_luminance_table_q50 = [16 11 10 16 24 40 51 61;
12 12 14 19 26 58 60 55;
14 13 16 24 40 57 69 56;
14 17 22 29 51 87 80 62;
18 22 37 56 68 109 103 77;
24 36 55 64 81 104 113 92;
49 64 78 87 103 121 120 101;
72 92 95 98 112 100 103 99;]

standard_chrominance_table_q50 = [17 18	24	47	99	99	99	99;
18	21	26	66	99	99	99	99;
24	26	56	99	99	99	99	99;
47	66	99	99	99	99	99	99;
99	99	99	99	99	99	99	99;
99	99	99	99	99	99	99	99;
99	99	99	99	99	99	99	99;
99	99	99	99	99	99	99	99;]

function JPEG_quality_multiplier(q)
    if q < 50
        mult = 50 / q
    else
        mult = (100 - q) / 50 
    end
    return round(mult, digits=2)
end

function get_luminance_quant_table(quality)
    return replace!(Int64.(round.(JPEG_quality_multiplier(quality) .* standard_luminance_table_q50)), 0 => 1)
end

function get_chrominance_quant_table(quality)
    return replace!(Int64.(round.(JPEG_quality_multiplier(quality) .* standard_chrominance_table_q50)), 0 => 1)
end

function find_nearest_quality_factor_luminance(mle_matrix, values_to_consider)
    nearest_dist = Inf
    nearest_q = -1
    for q in 1:99 #Just ignore the 100 case because that is the all 1's 
        quant_table = get_luminance_quant_table(q)
        mae = sum(abs.(mle_matrix[values_to_consider] .- quant_table[values_to_consider]))
        if mae < nearest_dist
            nearest_dist = mae
            nearest_q = q
        end
    end
    return nearest_q
end

function find_nearest_quality_factor_chrominance(mle_matrix, values_to_consider)
    nearest_dist = Inf
    nearest_q = -1
    for q in 1:99 #Just ignore the 100 case because that is the all 1's 
        quant_table = get_chrominance_quant_table(q)
        mae = sum(abs.(mle_matrix[values_to_consider] .- quant_table[values_to_consider]))
        if mae < nearest_dist
            nearest_dist = mae
            nearest_q = q
        end
    end
    return nearest_q
end

function DCT_coefficients(base_image)
    dct_coefficients = Int64.(zeros(length(base_image) ÷ 64, 8, 8))
    block_index = 1
    for i in 1:8:axes(base_image)[1][end]
        for j in 1:8:axes(base_image)[2][end]
            block = Float64.(base_image[i:i+7, j:j+7])
            block_coefficients = FFTW.dct(block .- 128)
            dct_coefficients[block_index, :, :] = Int64.(round.(block_coefficients))
            block_index += 1
        end
    end
    return dct_coefficients
end