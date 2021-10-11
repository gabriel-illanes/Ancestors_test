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
  recomb_offsrping::Vector{Float64} = [];
  pops_offsrping::Vector{String} = [];
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
      push!(hap_lenghts, recomb_sites[ind + j + 1] - recomb_sites[ind]);
      push!(hap_pops, pops[ind]);
      ind = ind + j + 1;
    end
  end
  return hap_lengths, hap_pops
end

function family_tree_admix!(recomb_tree::Vector{Vector{Float64}}, tree_pop::Vector{Vector{String}})
  len = length(recomb_tree);
  if !(largo == 1)
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
  pop_tree = Array{Array{String, 1}, 1}(undef, 2^t);
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
  hap_lenghts, hap_pops = get_hap_lenghts(recombs, pops);
  return sum(hap_pops .== target_pop) == 0 ? 0.0 : maximum(hap_lenghts[hap_pops .== target_pop])
end

function simulate_cdf_max_hap(max_iter::Int64, l::Float64, t::Int64, t_ancestors::Vector{Int64}, pₑ::Vector{Float64}, pₐ::Vector{Float64}, pₙ::Vector{Float64})
  sim = Vector{Float64}(undef,max_iter);
  for ind = 1:max_iter
    sim[ind] = max_hap(family_tree(l, t, t_ancestros, pₑ, pₐ, pₙ)...)
  end
  return sim
end

function sum_hap(recombs::Vector{Float64}, pops::Vector{String}, target_pop::String = "native")
  hap_lenghts, hap_pops = get_hap_lenghts(recombs, pops);
  return sum(hap_pops .== target_pop) == 0 ? 0.0 : sum(hap_lenghts[hap_pops .== target_pop])
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
  simu = rand(22, N)
  for i = 1:22
    for j = 1:N
      simu[i,j] = dist_pᵢ(simu[i,j], atoms_all[i, :])
    end
  end
  sim_pi = vec(sum(simu, dims = 1))
  return sum(sim_pi .< pm)/N
end


function test_maxmax(m_obs::Vector{Float64}, max_iter::Int64, t::Int64,
   L::Vector{Float64})::Float64
  t_ancestors = fill(0, 2^t);
  pₑ = push!(fill(1.0, 2^t - 1), 0.0);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0.0, 2^t - 1), 1.0);
  sim_all = simulate_max_all(max_iter, L, t, t_ancestros, pₑ, pₐ, pₙ)
  atomos_all = matriz_atomos(sim_all, L)
  pi_all = vector_pi(m_obs, sim_all, atomos_all)
  pvalue = calcular_pmax(pi_all, atomos_all)
  return pvalue
end

function test_maxsum(m_obs::Vector{Float64}, max_iter::Int64, t::Int64,
   L::Vector{Float64})::Float64
  t_ancestros = fill(0, 2^t);
  pₑ = push!(fill(1.0, 2^t - 1), 0.0);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0.0, 2^t - 1), 1.0);
  sim_all = simulaciones_max(max_iter, L, t, t_ancestors, pₑ, pₐ, pₙ)
  atomos_all = atoms_matrix(sim_all, L)
  pi_all = vector_pᵢ(m_obs, sim_all, atoms_all)
  pvalue = calcular_psum(pi_all, atomos_all)
  return pvalue
end


function test_summax(m_obs::Vector{Float64}, max_iter::Int64, t::Int64,
   L::Vector{Float64})::Float64
  t_ancestros = fill(0, 2^t);
  pₑ = push!(fill(1.0, 2^t - 1), 0.0);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0.0, 2^t - 1), 1.0);
  sim_all = simulaciones_sum(max_iter, L, t, t_ancestros, pₑ, pₐ, pₙ)
  atomos_all = matriz_atomos(sim_all, L)
  pi_all = vector_pi(m_obs, sim_all, atomos_all)
  pvalue = calcular_pmax(pi_all, atomos_all)
  return pvalue
end

function test_sumsum(m_obs::Vector{Float64}, max_iter::Int64, t::Int64,
   L::Vector{Float64})::Float64
  t_ancestros = fill(0, 2^t);
  pₑ = push!(fill(1.0, 2^t - 1), 0.0);
  pₐ = fill(0.0, 2^t);
  pₙ = push!(fill(0.0, 2^t - 1), 1.0);
  sim_all = simulaciones_sum(max_iter, L, t, t_ancestros, pₑ, pₐ, pₙ)
  atomos_all = matriz_atomos(sim_all, L)
  pi_all = vector_pi(m_obs, sim_all, atoms_all)
  pvalue = get_psum(pi_all, atoms_all)
  return pvalue
end

function test_method_all(M_obs::Matrix{Float64}, max_iter::Int64, T::Vector{Int64},
  L::Vector{Float64}, method::Function)::Matrix{Float64}
  N_donors = size(M_obs)[2]
  N_times = length(T)
  Pvalues = Matrix{Float64}(undef, N_donors, N_times)
  for t = 1:N_times
    for d = 1:N_donors
      Pvalues[d, t] = method(M_obs[:, d], max_iter, T[t], L)
    end
    print("t = ", T[t] + 1, " listo \n")
  end
  return Pvalues
end
