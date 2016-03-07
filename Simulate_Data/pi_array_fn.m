
function Pi_array = pi_array_fn( Pi_star_array, Zetaj_shocks_array, Zetajft_shocks_array, sigma_zetaj, sigma_zetajft)


Pi_array = Pi_star_array + sigma_zetaj * Zetaj_shocks_array + sigma_zetajft * Zetajft_shocks_array;

end
