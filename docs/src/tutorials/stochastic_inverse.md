# Stochastic inversion

```@setup rp_si
using Porosity, Distributions, Turing, CairoMakie
```

## Single parameter inference

Now that we know how to compute the rock physics responses, we move towards performing stochastic inversion for rock physics parameters given some rock conductivity and the corresponding uncertainty.

!!! note
    
    When performing stochastic inversion, make sure you output the vectors from the models, even if its a `(1,1)` or `(1,)` vector

```@example rp_si
# Creating synthetic data
m = Poe2010([1000.0], 100.0)
resp = forward(m, [])
err_resp = RockphyCond(0.01 .* abs.(resp.σ))
```

Then create an *a priori* using corresponding distribution, here `Poe2010Distribution`. In the above, we want to infer only for temperature and water content in olivine, so we'll define the distribution as such.

!!! note
    
    Note the resemblance between (`Poe2010`)[@ref] and (`Poe2010Distribution`)[@ref]

```@example rp_si
mdist = Poe2010Distribution(MvNormal([1400.0], [400.0;]), [100.0])
```

Then defining the likelihood

```@example rp_si
rdist = RockphyCondDistribution(Porosity.normal_dist)
```

and putting everything together to perform the inference using `stochastic inverse`[@ref]

```@example rp_si
m_cache = mcmc_cache(mdist, rdist, 1000, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 1000, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example rp_si
f = Figure()
ax = Axis(f[1, 1]; title="T(K)")

i = mcmc_chain_prior.name_map[1][1]

dat_prior = mcmc_chain_prior[i][:]
dat_posterior = mcmc_chain_posterior[i][:]

density!(ax, dat_prior; label="prior")
density!(ax, dat_posterior; label="posterior")

f[1, 2] = Legend(f, ax)
nothing # hide
```

```@raw html
</details>
```

```@example rp_si
f # hide
```

## Double parameter inference

Lets try to also infer the water content along with the temperature. Everything remains the same, except that now, in `ps_nt`, the vector corresponding to water content `Ch2o_ol` will be replaced by a corresponding distribution.

```@example rp_si
mdist = Poe2010Distribution(
    MvNormal([1200.0], [400.0;]), product_distribution([Uniform(50.0, 150.0)]))
rdist = RockphyCondDistribution(Porosity.normal_dist)

m_cache = mcmc_cache(mdist, rdist, 1000, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 1000, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example rp_si
f = Figure()
ax1 = Axis(f[1, 1]; title="T(K)")

i = mcmc_chain_prior.name_map[1][1]

dat_prior = mcmc_chain_prior[i][:]
dat_posterior = mcmc_chain_posterior[i][:]

density!(ax1, dat_prior; label="prior")
density!(ax1, dat_posterior; label="posterior")

ax2 = Axis(f[1, 2]; title="water conc. (ppm)")

i = mcmc_chain_prior.name_map[1][2]

dat_prior = mcmc_chain_prior[i][:]
dat_posterior = mcmc_chain_posterior[i][:]

density!(ax2, dat_prior; label="prior")
density!(ax2, dat_posterior; label="posterior")

f[2, 1:2] = Legend(f, ax2; orientation=:horizontal)
nothing # hide
```

```@raw html
</details>
```

```@example rp_si
f # hide
```

## Two-phase inference

Lets combine `SEO3` with `Sifre2014` and also try to infer the porosity.

```@example rp_si
# Creating synthetic data
m = two_phase_modelType(SEO3, Sifre2014, HS1962_plus)
ps_nt = (; T=[1400.0] .+ 273, Ch2o_m=[100.0], Cco2_m=[100.0], ϕ=[0.1])
model = m(ps_nt)
resp = forward(model, [])
err_resp = RockphyCond(0.01 .* abs.(resp.σ))

m = two_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution, HS1962_plus)
ps_nt_dist = (;
    T=MvNormal([1000.0], [400.0;]),
    Ch2o_m=product_distribution([Uniform(50.0, 150.0)]),
    Cco2_m=[100.0], ϕ=product_distribution([Uniform(0.01, 0.2)])
)
mdist = m(ps_nt_dist)
rdist = RockphyCondDistribution(Porosity.normal_dist)
```

```@example rp_si
m_cache = mcmc_cache(mdist, rdist, 1000, Prior());
mcmc_chain_prior = stochastic_inverse(resp, err_resp, [], m_cache);

m_cache = mcmc_cache(mdist, rdist, 1000, NUTS());
mcmc_chain_posterior = stochastic_inverse(resp, err_resp, [], m_cache);
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example rp_si
f = Figure()

labels = ["ϕ", "T (K)", "water conc. (ppm)"]

for i in 1:3
    axi = Axis(f[1, i]; title=labels[i])

    idx = mcmc_chain_prior.name_map[1][i]

    dat_prior = mcmc_chain_prior[idx][:]
    dat_posterior = mcmc_chain_posterior[idx][:]

    density!(axi, dat_prior; label="prior")
    density!(axi, dat_posterior; label="posterior")
end

f[2, 1:3] = Legend(f, f.content[end]; orientation=:horizontal)
nothing # hide
```

```@raw html
</details>
```

```@example rp_si
f # hide
```

## Multi rock physics inference

Things work similarly here, but we provide an example for convenience:

```@example rp_si
m = multi_rp_modelType(SEO3, anharmonic_poro, Nothing, Nothing)
ps_nt = (; T=[1000.0] .+ 273, P=3.0, ρ=3300.0, Ch2o_m=1000.0, ϕ=0.1)
model = m(ps_nt)
resp = forward(model, [])
err_resp = multi_rp_response(
    RockphyCond(0.01 .* abs.(resp.cond.σ)),
    RockphyElastic(
        0.01 .* resp.elastic.K,
        0.01 .* resp.elastic.G,
        0.01 .* resp.elastic.Vp,
        0.01 .* resp.elastic.Vs
    ),
    Nothing, Nothing
)

ps_nt_dist = (;
    T=product_distribution([Uniform(1200.0, 1400.0)]),
    P=[3.0],
    ρ=[3300.0],
    Ch2o_m=MvNormal([100.0], [20.0;]),
    ϕ=product_distribution([Uniform(0.01, 0.2)])
)

m1 = multi_rp_modelDistributionType(
    SEO3Distribution, anharmonic_poroDistribution, Nothing, Nothing)
mdist = m1(ps_nt_dist)
# rdist = RockphyCondDistribution(Porosity.normal_dist)
rdist = Porosity.multi_rp_responseDistribution(
    RockphyCondDistribution(normal_dist),
    RockphyElasticDistribution(normal_dist, normal_dist, normal_dist, normal_dist),
    Nothing, Nothing
)
nothing # hide
```

```@example rp_si
m_cache = mcmc_cache(mdist, rdist, 1000, Prior());
mcmc_chain_prior = stochastic_inverse(
    resp, err_resp, [], m_cache; response_fields=[:Vp, :Vs, :σ]);

m_cache = mcmc_cache(mdist, rdist, 1000, NUTS());
mcmc_chain_posterior = stochastic_inverse(
    resp, err_resp, [], m_cache; response_fields=[:Vp, :Vs, :σ]);
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example rp_si
f = Figure()

labels = ["T (K)", "ϕ"]

for i in 1:2
    axi = Axis(f[1, i]; title=labels[i])

    idx = mcmc_chain_prior.name_map[1][i]

    dat_prior = mcmc_chain_prior[idx][:]
    dat_posterior = mcmc_chain_posterior[idx][:]

    density!(axi, dat_prior; label="prior")
    density!(axi, dat_posterior; label="posterior")
end

f[2, 1:2] = Legend(f, f.content[end]; orientation=:horizontal)
nothing # hide
```

```@raw html
</details>
```

```@example rp_si
f # hide
```

### Multi rock physics with two phase

As a last example, lets take the two phases case above and see if we can do a better job at inference if we also had the `Vp`, `Vs` data

```@example rp_si
# Synthetic data
m_mix = two_phase_modelType(SEO3, Sifre2014, HS1962_plus)
m = multi_rp_modelType(typeof(m_mix), anharmonic, Nothing, Nothing)

ps_nt = (; T=[1200.0] .+ 273, Ch2o_m=[100.0], Cco2_m=[100.0], ϕ=[0.1], P=[3.0], ρ=[3300.0])

model = m(ps_nt)
resp = forward(model, [])

err_resp = multi_rp_response(
    RockphyCond(0.01 .* abs.(resp.cond.σ)),
    RockphyElastic(
        0.01 .* resp.elastic.K,
        0.01 .* resp.elastic.G,
        0.01 .* resp.elastic.Vp,
        0.01 .* resp.elastic.Vs
    ),
    Nothing, Nothing
)

# Defining apriori

ps_nt_dist = (;
    T=product_distribution([Uniform(1270.0, 1670.0)]),
    Ch2o_m=product_distribution([Uniform(50.0, 150.0)]),
    Cco2_m=[100.0], ϕ=product_distribution([Uniform(0.01, 0.2)]),
    P=[3.0],
    ρ=[3300.0]
)

m_mixdist = two_phase_modelDistributionType(
    SEO3Distribution, Sifre2014Distribution, HS1962_plus)
m = multi_rp_modelDistributionType(
    typeof(m_mixdist), anharmonicDistribution, Nothing, Nothing)

mdist = m(ps_nt_dist);
rdist = Porosity.multi_rp_responseDistribution(
    RockphyCondDistribution(normal_dist),
    RockphyElasticDistribution(normal_dist, normal_dist, normal_dist, normal_dist),
    Nothing, Nothing
)

m_cache = mcmc_cache(mdist, rdist, 1000, Prior());
mcmc_chain_prior = stochastic_inverse(
    resp, err_resp, [], m_cache; response_fields=[:Vp, :Vs, :σ]);

m_cache = mcmc_cache(mdist, rdist, 1000, NUTS());
mcmc_chain_posterior = stochastic_inverse(
    resp, err_resp, [], m_cache; response_fields=[:Vp, :Vs, :σ]);
nothing # hide
```

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example rp_si
f = Figure()

labels = ["ϕ", "T (K)", "water conc. (ppm)"]

for i in 1:3
    axi = Axis(f[1, i]; title=labels[i])

    idx = mcmc_chain_prior.name_map[1][i]

    dat_prior = mcmc_chain_prior[idx][:]
    dat_posterior = mcmc_chain_posterior[idx][:]

    density!(axi, dat_prior; label="prior")
    density!(axi, dat_posterior; label="posterior")
end

f[2, 1:3] = Legend(f, f.content[end]; orientation=:horizontal)
nothing # hide
```

```@raw html
</details>
```

```@example rp_si
f # hide
```
