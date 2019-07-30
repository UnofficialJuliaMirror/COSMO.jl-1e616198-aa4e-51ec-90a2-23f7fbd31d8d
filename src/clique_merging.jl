function print_merge_logs(ws; print_max = 10)
  sp_arr = ws.ci.sp_arr
  # print the merge logs for each sparsity pattern
  println(">>> Merge Logs:")
  for (iii, sp) in enumerate(sp_arr)
    m_log = sp.sntree.merge_log
    println("Sparsity Pattern Nr. $(iii), Graph Size: $(length(sp.sntree.par))")
    println("\t Num merges: $(m_log.num)\n\t Num decisions: $(length(m_log.decisions))\n\t Pairs:")

    # only print a certain number
    print_max = min(print_max, length(m_log.decisions))
    for jjj = 1:print_max
      p = m_log.clique_pairs[jjj, :]
      println("\t \tC_A: $(p[1]), C_B: $(p[2]) ($(m_log.decisions[jjj]))")
    end

  end
end



get_merge_logs(ws) = map(x -> x.sntree.merge_log, ws.ci.sp_arr)


# -------------------------------------
# Functions related to clique merging
# -------------------------------------
"""
    AbstractMergeStrategy

A merge strategy determines how the cliques in a clique tree are merged to improve computation time of the projection. Each merge receipt should
implement the following functions:

 - initialise!: Initialise the graph / tree from the input clique tree, allocate memory
 - traverse: A method that determines how the clique tree / graph is traversed
 - evaluate: A method that decides whether to merge two cliques
 - update!: A method to update local strategy-related information after a merge
"""
abstract type AbstractMergeStrategy end

"""
    NoMerge <: AbstractMergeStrategy

A strategy that does not merge cliques.
"""
 struct NoMerge <: AbstractMergeStrategy end

abstract type AbstractEdgeScore end


struct RelIntersect <: AbstractEdgeScore end
struct ComplexityScore <: AbstractEdgeScore end
merge_cliques!(t, strategy) = merge_cliques!(t, strategy())
merge_cliques!(t::SuperNodeTree, strategy::NoMerge) = nothing
function merge_cliques!(t::SuperNodeTree, strategy::AbstractMergeStrategy)
  # strategy = strategy()

  initialise!(t, strategy)

  while !strategy.stop
    # find two merge candidates
    cand = traverse(t, strategy)
    ordered_cand = copy(cand)
    # evaluate wether to merge the candidates
    do_merge = evaluate(t, strategy, cand)
    if do_merge
      # return cand in such a way that the remaining clique comes first
      ordered_cand[:, :] = merge_two_cliques!(t, cand, strategy)
      # decrement number of cliques in tree
      t.num -= 1
    end
    log_merge!(t, do_merge, ordered_cand)

    t.num == 1 && break
    strategy.stop && break
    # update strategy information after the merge
    update!(strategy, t, cand, ordered_cand)
  end

  # number of cliques
  Nc = t.num
  # the merging creates empty supernodes and seperators, recalculate a post order for the supernodes
  snd_post = post_order(t.snd_par, t.snd_child, Nc)
  t.snd_post = snd_post

  return nothing
end

# calculate block sizes (notice: the block sizes are stored in post order)
function calculate_block_dimensions!(t::SuperNodeTree)
  Nc = t.num
  t.nBlk = zeros(Nc)
  for iii = 1:Nc
    c = t.snd_post[iii]
    t.nBlk[iii] = Base.power_by_squaring(length(t.sep[c]) + length(t.snd[c]), 2)
  end
end

function log_merge!(t::SuperNodeTree, do_merge::Bool, cand::Array{Int64, 1})
  t.merge_log.clique_pairs = vcat(t.merge_log.clique_pairs, cand')
  push!(t.merge_log.decisions, do_merge)
  do_merge && (t.merge_log.num += 1)
  return nothing
end




# Merge two cliques that are in a parent - child relationship
function merge_child!(t::SuperNodeTree, cand::Array{Int64, 1})
  # determine which clique is the parent
  p, ch = determine_parent(t, cand[1], cand[2])

  # merge child's vertex sets into parent's vertex set
  push!(t.snd[p], t.snd[ch]...)
  t.snd[ch] = [] #necessary or just leave it
  t.sep[ch] = []

  # update parent structure
  @. t.snd_par[t.snd_child[ch]] = p
  t.snd_par[ch] = -1 #-1 instead of NaN, effectively remove that entry from the parent list

  # update children structure
  filter!(x -> x != ch, t.snd_child[p])
  push!(t.snd_child[p], t.snd_child[ch]...)
  t.snd_child[ch] = []
  return [p; ch]
end

# Merge two cliques that are in a sibling relationship
function merge_sibling!(t::SuperNodeTree, cand::Array{Int64, 1})
  c1 = cand[1]
  c2 = cand[2]
  # merge vertex set of cand[2] into cand[1]
  push!(t.snd[c1], t.snd[c2]...)
  t.snd[c2] = []
  union!(t.sep[c1], t.sep[c2])
  t.sep[c2] = []

  @. t.snd_par[t.snd_child[c2]] = c1
  t.snd_par[c2] = -1

  # update children structure
  push!(t.snd_child[c1], t.snd_child[c2]...)
  t.snd_child[c2] = []

  return [c1; c2]
end




"""
    PairwiseMerge <: AbstractMergeStrategy

A merge strategy that calculates the edge metric `A ∩ B / A ∪ B` for every two cliques that are in a parent-child or sibling relationship. The resulting clique
graph is traversed from the highest edge metric to the lowest.
"""
mutable struct PairwiseMerge <: AbstractMergeStrategy
  stop::Bool
  edges::AbstractMatrix
  edge_score::AbstractEdgeScore
  function PairwiseMerge(; edge_score = RelIntersect())
    new(false, spzeros(Float64, 0, 0), edge_score)
  end
end

# compute the edge set of the initial clique graph, only consider parent-child and sibling-sibling relationships
function initialise!(t, strategy::PairwiseMerge)
  n_cliques = length(t.snd_par)
  strategy.edges = spzeros(Float64, n_cliques, n_cliques)

  # following a descending post-order, add edges
  for i = length(t.snd_post):-1:1
    c = t.snd_post[i]
    children = t.snd_child[c]
    n_children = length(children)
    # add edge to each child  (Upper Triangle)
    for (j, ch) in enumerate(children)
      strategy.edges[min(c, ch),max(c, ch)] = edge_metric_parent(t.snd, t.sep, c, ch, strategy.edge_score)

      # add edges between siblings combinations (Lower Triangle)
      for sib in view(children,(j+1):n_children)
        strategy.edges[max(sib, ch), min(c, ch)] = edge_metric_siblings(t.snd, t.sep, ch, sib, strategy.edge_score)
      end
    end
  end
  return nothing
end

function traverse(t, strategy::PairwiseMerge)
  # find maximum edge value in sparse edge matrix
  return max_elem(strategy.edges)
end

# Find the matrix indices (i, j) of the first maximum element among the elements stored in A.nzval
function max_elem(A::SparseMatrixCSC)
  length(A.nzval) == 0 && throw(DomainError("Sparse matrix A doesn't contain any entries"))
  ~, ind = findmax(A.nzval)
  r = A.rowval[ind]

  c = findfirst(x -> x >= ind, A.colptr)
  # zero columns appear with repetitive indices in colptr, move behind the last zero column
  if c < length(A.colptr)
    while A.colptr[c] == A.colptr[c + 1]
     # global c;
      c += 1
    end
  end
  # col_ptr: 1 2 5
  # if column 2 has range 2:5 and ind = 4, then c would be returned as 3, but you want the start, i.e. 2
  if A.colptr[c] > ind
    c -= 1
  end
  return [r; c]
end


merged_block_size(t::SuperNodeTree, c1::Int64, c2::Int64) = c1 < c2 ? block_size_child(t, c1, c2) : block_size_sibling(t, c1, c2)

# Given two cliques c1 and c2, return the parent clique first
function determine_parent(t::SuperNodeTree, c1::Int64, c2::Int64)
  if in(c2, t.snd_child[c1])
    return c1, c2
  else
    return c2, c1
  end
end

function block_size_child(t::SuperNodeTree, c1::Int64, c2::Int64)
  p, ch = determine_parent(t, c1, c2)
  return length(t.snd[p]) + length(t.sep[p]) + length(t.snd[ch])
end

function block_size_sibling(t::SuperNodeTree, c1::Int64, c2::Int64)
  return length(t.snd[c1]) + length(t.snd[c2]) + union_dim(t.sep[c1], t.sep[c2])
end

evaluate(t, strategy::PairwiseMerge, cand) = evaluate(t, strategy, cand, strategy.edge_score)


function evaluate(t, strategy::PairwiseMerge, cand, edge_score::RelIntersect)

  c1 = cand[1]
  c2 = cand[2]

  # case if only one clique remains
  if c1 == c2
    strategy.stop = true
    return false
  end

  # individual block sizes
  n_1 = length(t.snd[c1]) + length(t.sep[c1])
  n_2 = length(t.snd[c2]) + length(t.sep[c2])


  n_m = merged_block_size(t, c1, c2)
  n_ops_diff = compute_complexity_savings(n_1, n_2, n_m)

  do_merge = (n_ops_diff >= 0)

  if !do_merge
    strategy.stop = true
  end
  return do_merge
end

function evaluate(t, strategy::PairwiseMerge, cand, edge_score::ComplexityScore)
  do_merge = (strategy.edges[cand[1], cand[2]] >= 0)

  if !do_merge
    strategy.stop = true
  end
  return do_merge
end

# Assuming the complexity of the projection is roughly O(n^3), how many operations are saved by projection the merged cliques
# instead of the individual cliques
compute_complexity_savings(n_1::Int64, n_2::Int64, n_m::Int64) = n_1^3 + n_2^3 - n_m^3

# Approximates the number of operations for one projection of all the cliques in the tree
compute_complexity(t::COSMO.SuperNodeTree) = sum(map(x -> x^3, t.nBlk))


function merge_two_cliques!(t::SuperNodeTree, cand::Array{Int64, 1}, strategy::PairwiseMerge)
  # parent - child relationships are stored in upper triangle
  cand[1] < cand[2] ? merge_child!(t, cand) : merge_sibling!(t, cand)
end

# After a merge operation update the information of the strategy
function update!(strategy::PairwiseMerge, t, cand, ordered_cand)
  c1 = ordered_cand[1]
  c_removed = ordered_cand[2]
  edges = strategy.edges

  # remove edge value between c1 and c2 (in order) ...
  edges[cand[1], cand[2]] = 0
  # ... and all to the removed clique
  edges[c_removed, :] .= 0
  edges[:, c_removed] .= 0

  # recalculate edge values of c1's children
  for ch in t.snd_child[c1]
    edges[min(c1, ch), max(c1, ch)] = COSMO.edge_metric_parent(t.snd, t.sep, c1, ch, strategy.edge_score)
  end

  # edge to c1's parent (for all non-root nodes)
  p = t.snd_par[c1]
  if p > 0
    edges[min(p, c1), max(p, c1)] = COSMO.edge_metric_parent(t.snd, t.sep, p, c1, strategy.edge_score)

    # edges to c1's siblings
    siblings = t.snd_child[p]
    for sib in siblings
      if sib != c1
          edges[max(c1, sib), min(c1, sib)] = COSMO.edge_metric_siblings(t.snd, t.sep, c1, sib, strategy.edge_score)
      end
    end
  end
  dropzeros!(strategy.edges)
end


# computes the edge metric for cliques C_a and C_b: (C_a ∩ C_b) / (C_a ∪ C_b) where C_a is the parent of C_b
function edge_metric_parent(res, sep, c_a, c_b, edge_score::RelIntersect)
  return length(sep[c_b]) / (length(sep[c_a]) + length(res[c_a]) + length(res[c_b]))
end

# computes the edge metric for cliques C_a and C_b: (C_a ∩ C_b) / (C_a ∪ C_b) where C_a is the sibling of C_b
function edge_metric_siblings(res, sep, c_a, c_b, edge_score::RelIntersect)
  return intersect_dim(sep[c_a], sep[c_b]) / (length(res[c_a]) + length(res[c_b]) + union_dim(sep[c_a], sep[c_b]))
end

# computes the edge metric for cliques C_a and C_b in terms of how the number of projection operationss change when merged
function edge_metric_parent(res, sep, c_a, c_b, edge_score::ComplexityScore)

  n_1 = length(res[c_a]) + length(sep[c_a])
  n_2 = length(res[c_b]) + length(sep[c_b])
  # merged block size
  n_m = n_1 + length(res[c_b])
  return compute_complexity_savings(n_1, n_2, n_m)
end

# computes the edge metric for cliques C_a and C_b: (C_a ∩ C_b) / (C_a ∪ C_b) where C_a is the sibling of C_b
function edge_metric_siblings(res, sep, c_a, c_b, edge_score::ComplexityScore)
  n_1 = length(res[c_a]) + length(sep[c_a])
  n_2 = length(res[c_b]) + length(sep[c_b])
  # merged block size
  n_m = length(res[c_a]) + length(res[c_b]) + union_dim(sep[c_a], sep[c_b])
  return compute_complexity_savings(n_1, n_2, n_m)
end

# Finds the size of the set A ∩ B under the assumption that B only has unique elements
function intersect_dim(A::AbstractVector, B::AbstractVector)
  dim = 0
  for elem in B
    in(elem, A) && (dim += 1)
  end
  return dim
end

# Finds the size of the set A ∪ B under the assumption that A and B only have unique elements
function union_dim(A::AbstractVector, B::AbstractVector)
  dim = length(A)
  for elem in B
    !in(elem, A) && (dim += 1)
  end
  return dim
end

"""
    TreeTraversalMerge(t_fill = 5, t_size = 5) <: AbstractMergeStrategy

The merge strategy suggested in Sun / Andersen - Decomposition in conic optimization with partially separable structure (2013).
The initial clique tree is traversed in topological order and cliques are greedily merged to their parent if evaluate() returns true.
"""
mutable struct TreeTraversalMerge <: AbstractMergeStrategy
  stop::Bool
  clique_ind::Int64
  t_fill::Int64
  t_size::Int64

  function TreeTraversalMerge(; t_fill = 5, t_size = 5)
    new(false, 2, t_fill, t_size)
  end
end

function initialise!(t, strategy::TreeTraversalMerge)
  # start with node that has second highest order
  strategy.clique_ind = length(t.snd) - 1
end
# traverse tree in descending topological order and return clique and its parent, root has highest order
function traverse(t, strategy::TreeTraversalMerge)
  c = t.snd_post[strategy.clique_ind]
  return [t.snd_par[c]; c]
end

function evaluate(t, strategy::TreeTraversalMerge, cand)
  strategy.stop && return false
  c = cand[1]
  par = cand[2]
  n_sep = length(t.sep[c])
  n_sep_par = length(t.sep[par])
  n_cl = length(t.snd[c]) + n_sep
  n_cl_par = length(t.snd[par]) + n_sep
  return (n_cl_par - n_sep) * (n_cl - n_sep) <= strategy.t_fill || max(n_cl - n_sep, n_cl_par - n_sep_par) <= strategy.t_size
end

# TreeTraversal merge strategy always merges parent-child pairs
merge_two_cliques!(t::SuperNodeTree, cand::Array{Int64, 1}, strategy::TreeTraversalMerge) = merge_child!(t, cand)


function update!(strategy::TreeTraversalMerge, t, cand, ordered_cand)
  # try to merge last node of order 1, then stop
  if strategy.clique_ind == 1
    strategy.stop = true
  # otherwise decrement node index
  else
    strategy.clique_ind -= 1
  end
end