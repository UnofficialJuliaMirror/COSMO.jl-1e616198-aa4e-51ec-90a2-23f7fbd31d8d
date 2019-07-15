# Unit tests for clique merging
using COSMO, SparseArrays, LinearAlgebra, Test


mutable struct SnodeTree
  res::Array{Array{Int64, 1},1}
  sep::Array{Array{Int64, 1},1}
  par::Array{Int64}
  post::Array{Int64}
  child::Array{Array{Int64, 1},1}
  function SnodeTree(res, sep, par, post)
    # compute children structure
    child = COSMO.child_from_par(par)
    new(res, sep, par, post, child)
  end
end

# lets consider the following graph with cliques
res1 = [[1; 2], [3], [4]];
sep1 = [[], [1], [2]];
par1 = [0; 1; 1];
post1 = [1; 2; 3]
ref_edges = spzeros(3, 3)
ref_edges[1, 2] = 1 / 3;
ref_edges[1, 3] = 1 / 3;
ref_edges[3, 2] = 0;


tree1 = SnodeTree(res1, sep1, par1, post1)

@testset "Clique merging" begin

  # compute and verify the initial (approximate) "clique graph" layout and edge values (only consider parent-child and sibling-sibling)
  edges = COSMO.initial_clique_graph(tree1)
  @test edges == ref_edges

  cand = COSMO.find_merge_candidates(edges)
  @test cand == [1; 2]

  COSMO.merge_cliques!(tree1, edges, cand)
  @test tree1.res == [[1; 2; 3], [4]]
  @test tree1.sep == [[], [1], [2]]
  @test tree1.par == [0, 0, 1]
  res_edges[2, 1] = 0
  res_edge[3, 1] = 1 / 4;
  @test edges == res_edges
end
