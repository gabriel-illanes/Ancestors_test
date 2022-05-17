using Distributions
using HypothesisTests
using Random
using Plots

function create_recomb_sites(l::Float64, t::Int64)
  recomb_sites::Vector{Float64} = [0.0];
  if !(t == 0)
    while (site = recomb_sites[end] + rand(Exponential(1/t))) < l
      push!(recomb_sites, site);
    end
  end
  push!(recomb_sites, l);
  return recomb_sites
end

function create_ancestor(l::Float64, t::Int64, pₑ::Float64, pₐ::Float64, pₙ::Float64)

  recomb_sites = create_recomb_sites(l, t);
  pops::Vector{String} = wsample(["european", "african", "native"], [pₑ, pₐ, pₙ], length(recomb_sites)-1);
  return recomb_sites, pops
end


function admix!(recomb_offspring::Vector{Float64}, pops_offspring::Vector{String}, recomb_ancestors::Array{Array{Float64,1},1}, pops_ancestors::Array{Array{String,1},1}, new_recomb_sites::Vector{Float64}, start₁::Int64, start₂::Int64)
  if length(new_recomb_sites) == 1
    append!(recomb_offspring, recomb_ancestors[1][start₁:end]);
    append!(pops_offspring, pops_ancestors[1][max((start₁-1),1):end]);
  else
    let ind₁ = searchsortedlast(recomb_ancestors[1], new_recomb_sites[1]), ind₂ = searchsortedfirst(recomb_ancestors[2], new_recomb_sites[1])
      append!(recomb_offspring, recomb_ancestors[1][start₁:ind₁]);
      push!(recomb_offspring, new_recomb_sites[1]);
      append!(pops_offspring, pops_ancestors[1][max((start₁-1),1):ind₁]);
      admix!(recomb_offspring, pops_offspring, recomb_ancestors[[2,1]], pops_ancestors[[2,1]], new_recomb_sites[2:end], ind₂, ind₁ + 1);
    end
  end
end

function create_admix(recomb_ancestors::Array{Array{Float64,1},1}, pops_ancestors::Array{Array{String,1},1})
  recomb_offspring::Vector{Float64} = [];
  pops_offspring::Vector{String} = [];
  let start₁ = 1, start₂ = 1, l = recomb_ancestors[1][end], flip_coin = shuffle(1:2)
    new_recomb_sites = create_recomb_sites(l, 1);
    admix!(recomb_offspring, pops_offspring, recomb_ancestors[flip_coin], pops_ancestors[flip_coin], new_recomb_sites[2:end], 1, 1);
  end
  return recomb_offspring, pops_offspring
end


function get_hap_lengths(recomb_sites::Vector{Float64}, pops::Vector{String})
  hap_lengths = Vector{Float64}(undef, 0);
  hap_pops = Vector{String}(undef, 0);
  let ind = 1, l = length(recomb_sites)
    while ind < l
      j = 0;
      while ((ind + j + 1) <= (l-1))  && (pops[ind+j+1] == pops[ind])
        j +=1;
      end
      push!(hap_lengths, recomb_sites[ind + j + 1] - recomb_sites[ind]);
      push!(hap_pops, pops[ind]);
      ind = ind + j + 1;
    end
  end
  return hap_lengths, hap_pops
end

function family_tree_admix!(recomb_tree::Vector{Vector{Float64}}, tree_pop::Vector{Vector{String}})
  len = length(recomb_tree);
  if !(len == 1)
    for ind = 1:div(len,2)
      let son = create_admix(collect((pop!(recomb_tree), pop!(recomb_tree))), collect((pop!(tree_pop), pop!(tree_pop))))
        pushfirst!(recomb_tree, son[1]);
        pushfirst!(tree_pop, son[2]);
      end
    end
    family_tree_admix!(recomb_tree, tree_pop);
  end
end

function family_tree(l::Float64, t::Int64, t_ancestors::Vector{Int64}, pₑ::Vector{Float64}, pₐ::Vector{Float64}, pₙ::Vector{Float64})
  recomb_tree = Array{Array{Float64, 1}, 1}(undef, 2^t);
  tree_pop = Array{Array{String, 1}, 1}(undef, 2^t);
  for ind = 1:2^t
    let ancestor = create_ancestor(l, t_ancestors[ind], pₑ[ind], pₐ[ind], pₙ[ind])
      recomb_tree[ind] = ancestor[1];
      tree_pop[ind] = ancestor[2];
    end
  end
  family_tree_admix!(recomb_tree, tree_pop);
  return pop!(recomb_tree), pop!(tree_pop)
end

function max_hap(recombs::Vector{Float64}, pops::Vector{String}, target_pop::String = "native")
  hap_lengths, hap_pops = get_hap_lengths(recombs, pops);
  return sum(hap_pops .== target_pop) == 0 ? 0.0 : maximum(hap_lengths[hap_pops .== target_pop])
end

function simulate_cdf_max_hap(max_iter::Int64, l::Float64, t::Int64, t_ancestors::Vector{Int64}, pₑ::Vector{Float64}, pₐ::Vector{Float64}, pₙ::Vector{Float64})
  sim = Vector{Float64}(undef,max_iter);
  for ind = 1:max_iter
    sim[ind] = max_hap(family_tree(l, t, t_ancestors, pₑ, pₐ, pₙ)...)
  end
  return sim
end

function sum_hap(recombs::Vector{Float64}, pops::Vector{String}, target_pop::String = "native")
  hap_lengths, hap_pops = get_hap_lengths(recombs, pops);
  return sum(hap_pops .== target_pop) == 0 ? 0.0 : sum(hap_lengths[hap_pops .== target_pop])
end

function simulate_cdf_sum_hap(max_iter::Int64, l::Float64, t::Int64, t_ancestors::Vector{Int64}, pₑ::Vector{Float64}, pₐ::Vector{Float64}, pₙ::Vector{Float64})
  sim = Vector{Float64}(undef, max_iter);
  for ind = 1:max_iter
    sim[ind] = sum_hap(family_tree(l, t, t_ancestors, pₑ, pₐ, pₙ)...);
  end
  return sim
end


####

function get_atoms(sim::Vector{Float64}, l::Float64)::Vector{Float64}
  N = length(sim)
  p0 = sum(sim .== 0.)/N
  pL = sum(sim .== l)/N
  return [p0, pL]
end

function get_p(sim::Vector{Float64}, mi::Float64)::Float64
  N = length(sim)
  return sum(sim .≤ mi)/N
end

function get_pᵢ(p::Float64, atoms::Vector{Float64})
  return (p-atoms[1])/(1.0 - atoms[1])
end

###
# Functions that simulate statistics under the null hypothesis,
# for both maximum and sum statistics
# p_e, p_a and p_n are vectors that represent average percentage of genetic
# information of each ancestral population

function simulate_max_all(max_iter::Int64, L::Vector{Float64}, t::Int64, t_ancestors::Vector{Int64},
        pₑ::Vector{Float64}, pₐ::Vector{Float64}, pₙ::Vector{Float64})::Matrix{Float64}
    sim_all = Matrix{Float64}(undef, 22, max_iter)
    for i = 1:22
      sim_all[i, :] = simulate_cdf_max_hap(max_iter, L[i], t, t_ancestors, pₑ, pₐ, pₙ)
    end
    return sim_all
end

function simulate_sum_all(max_iter::Int64, L::Vector{Float64}, t::Int64, t_ancestors::Vector{Int64},
        pₑ::Vector{Float64}, pₐ::Vector{Float64}, pₙ::Vector{Float64})::Matrix{Float64}
    sim_all = Matrix{Float64}(undef, 22, max_iter)
    for i = 1:22
      sim_all[i, :] = simulate_cdf_sum_hap(max_iter, L[i], t, t_ancestors, pₑ, pₐ, pₙ)
    end
    return sim_all
end

###

function atoms_matrix(sim_all::Matrix{Float64}, L::Vector{Float64})
   atoms_all = Matrix{Float64}(undef, 22, 2)
   for i = 1:22
     atoms_all[i, :] = get_atoms(sim_all[i, :], L[i])
   end
   return atoms_all
 end

function vector_pᵢ(m_obs::Vector{Float64}, sim_all::Matrix{Float64},
  atoms_all::Matrix{Float64})::Vector{Float64}
  p_all = Vector{Float64}(undef, 22)
  for i = 1:22
    p_all[i] = get_p(sim_all[i, :], m_obs[i])
  end
  pi_all = (p_all - atoms_all[:, 1])./(1.0 .- atoms_all[:, 1])
  return pi_all
end

###
# Functions to obtain the test p-value from all chromosome scores

function get_pmax(pi_all::Vector{Float64}, atoms_all::Matrix{Float64})::Float64
  pm = maximum(pi_all)
  pvalue = 1.
  for i = 1:22
    pvalue *= atoms_all[i,1] + pm * (1.0 - atoms_all[i,1] - atoms_all[i,2])
  end
  return pvalue
end

function get_psum(pi_all::Vector{Float64}, atoms_all::Matrix{Float64})::Float64
  function dist_pᵢ(u::Float64, atoms::Vector{Float64})::Float64
    if u < 1. - atoms[1] - atoms[2]
      return(u / (1. - atoms[1] - atoms[2]))
    elseif u < 1. - atoms[2]
      return 0.
    else
      return 1.
    end
  end

  N = 10000
  pm = sum(pi_all)
  sim = rand(22, N)
  for i = 1:22
    for j = 1:N
      sim[i,j] = dist_pᵢ(sim[i,j], atoms_all[i, :])
    end
  end
  sim_pi = vec(sum(sim, dims = 1))
  return sum(sim_pi .< pm)/N
end

###

function test_maxmax(m_obs::Vector{Float64}, max_iter::Int64, t::Int64,
   L::Vector{Float64})::Float64
  t_ancestors = fill(0, 2^t);
  pₑ = push!(fill(1.0, 2^t - 1), 0.0);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0.0, 2^t - 1), 1.0);
  sim_all = simulate_max_all(max_iter, L, t, t_ancestors, pₑ, pₐ, pₙ)
  atoms_all = atoms_matrix(sim_all, L)
  pi_all = vector_pᵢ(m_obs, sim_all, atoms_all)
  pvalue = get_pmax(pi_all, atoms_all)
  return pvalue
end

function test_maxsum(m_obs::Vector{Float64}, max_iter::Int64, t::Int64,
   L::Vector{Float64})::Float64
  t_ancestors = fill(0, 2^t);
  pₑ = push!(fill(1.0, 2^t - 1), 0.0);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0.0, 2^t - 1), 1.0);
  sim_all = simulate_max_all(max_iter, L, t, t_ancestors, pₑ, pₐ, pₙ)
  atoms_all = atoms_matrix(sim_all, L)
  pi_all = vector_pᵢ(m_obs, sim_all, atoms_all)
  pvalue = get_psum(pi_all, atoms_all)
  return pvalue
end


function test_summax(m_obs::Vector{Float64}, max_iter::Int64, t::Int64,
   L::Vector{Float64})::Float64
  t_ancestors = fill(0, 2^t);
  pₑ = push!(fill(1.0, 2^t - 1), 0.0);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0.0, 2^t - 1), 1.0);
  sim_all = simulate_sum_all(max_iter, L, t, t_ancestors, pₑ, pₐ, pₙ)
  atoms_all = atoms_matrix(sim_all, L)
  pi_all = vector_pᵢ(m_obs, sim_all, atoms_all)
  pvalue = get_pmax(pi_all, atoms_all)
  return pvalue
end

function test_sumsum(m_obs::Vector{Float64}, max_iter::Int64, t::Int64,
   L::Vector{Float64})::Float64
  t_ancestors = fill(0, 2^t);
  pₑ = push!(fill(1.0, 2^t - 1), 0.0);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0.0, 2^t - 1), 1.0);
  sim_all = simulate_sum_all(max_iter, L, t, t_ancestors, pₑ, pₐ, pₙ)
  atoms_all = atoms_matrix(sim_all, L)
  pi_all = vector_pᵢ(m_obs, sim_all, atoms_all)
  pvalue = get_psum(pi_all, atoms_all)
  return pvalue
end

##### test_method_all is the main function to run any test
# This function depends on all previous functions.
# method can take values in: test_maxmax, test_maxsum, test_summax, test_sumsum,
# depending on the choices of chromosome statistic and scores grouping
function test_method_all(M_obs::Matrix{Float64}, max_iter::Int64, T::Vector{Int64},
  L::Vector{Float64}, method::Function)::Matrix{Float64}
  N_donors = size(M_obs)[2]
  N_times = length(T)
  Pvalues = Matrix{Float64}(undef, N_donors, N_times)
  for t = 1:N_times
    for d = 1:N_donors
      Pvalues[d, t] = method(M_obs[:, d], max_iter, T[t], L)
    end
    print("t = ", T[t] + 1, " done \n")
  end
  return Pvalues
end

#####
## Importing the data
#####

using  DelimitedFiles
cd("/home/illa/Nextcloud/URUGENOMES/illa/Pres1")

# est_max has maximum of all lengths for each chromosome
est_max = readdlm("est_max.csv", ' ', Float64);

# Now we take the maximum among pairs of chromosomes,
# we obtain a statistic for each chromosome pair
est_max_fin = Matrix{Float64}(undef, 22, 20);
for i = 1:20
  est_max_fin[:, i] = max.(est_max[2*i-1, :], est_max[2*i, :])
end

# Now, the same thing for the sum of all lengths for each chromosome
est_sum = readdlm("est_sum.csv", ' ', Float64);
est_sum_fin = Matrix{Float64}(undef, 22, 20);
for i = 1:20
  est_sum_fin[:, i] = max.(est_sum[2*i-1, :], est_sum[2*i, :])
end

# L will have the length of every chromosome pair (in centimorgans)
# We scale to work in morgans
L = vec(readdlm("largos_chr.csv", ' ', Float64)./100);

#####
# Plots for test p-values for every individual

res = test_method_all(est_max_fin, 10000, [1, 2,3, 4], L, test_maxmax)
scatter(res, ylim = [0,1], label = ["t = 2" "t = 3" "t = 4" "t = 5"], legend = :topleft,
  xlabel = "Individuals", ylabel = "p-value", mode = "markers", marker = [:diamond :star :circle :square])
plot!([1,20], [0.05, 0.05], width = 2, label = "p-value = 0.05")

#####
# Power of the test: first scenario

# t is number of generations minus one
t = 4;
# We set the structure of anccestors.
t_ancestors = fill(4, 2^t);

# we set the chromosome coloring probabilities
pₑ = fill(1.0 - 2.0/2^t, 2^t);
pₐ = fill(0.0, 2^t);
pₙ = fill(2.0/2^t, 2^t);

# We simulate observations under this scenario
n_obs = 100
sims_sum = Matrix{Float64}(undef, 22, n_obs)
for i = 1:22
  sims_sum[i,:] = simulate_cdf_sum_hap(n_obs, L[i], t, t_ancestors, pₑ, pₐ, pₙ);
end

# We compute p-values for the chosen methodology
max_iter = 1000;
res = test_method_all(sims_sum, max_iter, t, L, test_summax)
sum(res .< 0.05, dims = 1)/max_iter


#####
# Power of the test: second scenario

# sample of chromosomes that will not have native american genetic information
chr_bool = [fill(1, 11); fill(0, 11)][randperm(22)]

max_iter = 1000
res_maxmax = Matrix{Float64}(undef, max_iter, 4);
for t = 1:4
  t_ancestors = fill(0, 2^t);
  sims_max = Matrix{Float64}(undef, 22, max_iter)
  for i = 1:22
    if chr_bool[i] < 0.5
      pₑ = push!(fill(1.0 , 2^t - 1), 0.);
      pₐ = fill(0.0, 2^t);
      pₙ = push!(fill(0., 2^t - 1), 1.);
    else
      pₑ = fill(1.0 , 2^t);
      pₐ = fill(0.0, 2^t);
      pₙ = fill(0., 2^t);
    end
    sims_max[i,:] = simulate_cdf_max_hap(max_iter, L[i], t, t_ancestors, pₑ, pₐ, pₙ);
  end
  res_maxmax[:, t] = test_method_all(sims_max, 1000, [t], L, test_maxmax)
end

sum(res_maxmax .< 0.05, dims = 1)/1000

#####
# Power of the test: third scenario

max_iter = 1000
res_summax = Matrix{Float64}(undef, max_iter, 4);
for t = 2:5
  t_ancestors = fill(0, 2^t);
  pₑ = push!(fill(1.0 , 2^t - 1), 0.);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0., 2^t - 1), 1.);

  sims_sum = Matrix{Float64}(undef, 22, max_iter)
  for i = 1:22
    sims_sum[i,:] = simulate_cdf_sum_hap(max_iter, L[i], t, t_ancestors, pₑ, pₐ, pₙ);
  end
  res_summax[:, t-1] = test_method_all(sims_sum, 1000, [t - 1], L, test_summax)
end

sum(res_summax .< 0.05, dims = 1)/1000
