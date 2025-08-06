fO2(T) = 10.0f0^(-24441.9f0 * inv(T) + 13.296f0)

arrh_dry(S, H, k, T) = S * exp(-H * inv(k * T))
arrh_wet(S, H, k, T, w, a, r) = S * (w^r) * exp(-(H - a * (w^inv(3.0f0))) * inv(k * T))
