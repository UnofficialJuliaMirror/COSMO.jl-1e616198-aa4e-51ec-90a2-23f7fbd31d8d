
mutable struct SuperNodeTree
	snd::Array{Array{Int64,1},1} #vertices of supernodes stored in one array (also called residuals)
	snd_par::Array{Int64,1}  # parent of supernode k is supernode j=snd_par[k]
	snd_post::Array{Int64,1} # post order of supernodal elimination tree
	snd_child::Array{Array{Int64,1},1}
	post::Array{Int64} # post ordering of the vertices in elim tree σ(j) = v
	par::Array{Int64}
	sep::Array{Array{Int64,1},1} #vertices of clique seperators
	nBlk::Array{Int64,1} #sizes of submatrizes defined by each clique, sorted by post-ordering, e.g. size of clique with order 3 => nBlk[3]
	function SuperNodeTree(L)
		par = etree(L)
		child = child_from_par(par)
		post = post_order(par, child)
		# sorting of children from highest order one to lowest make sure that a node always gets
		# added to the supernode with the highest rep vertex
		sort_children!(child, post)
		# faster algorithm according to Vandenberghe p.308
		degrees = higher_degrees(L)
		snd, snd_par = find_supernodes(par, child, post, degrees)

		snd_child = child_from_par(snd_par)
		# # a post-ordering of the elimination tree is needed to make sure that in the
		# # psd completion step the matrix is completed from lower-right to top-left
	 	snd_post = post_order(snd_par, snd_child)
	 	# given the supernodes (clique residuals) find the clique separators
		sep = find_separators(L, snd, snd_par, post)

		new(snd, snd_par, snd_post, snd_child, post, par, sep, [0])
	end
	# FIXME: only for debugging purposes
	function SuperNodeTree(res, par, post, sep)
		child = child_from_par(par)
  	new(res, par, post, child, [1], [1], sep, [1])
	end
end

function sort_children!(child, post)
	for v = 1:length(child)
		child[v] = sort(child[v], by = x -> post[x], rev = true)
	end
end
# -------------------------------------
# FUNCTION DEFINITIONS
# -------------------------------------
# given v=σ^-1(i) it returns i=σ(v)
function invert_order(sigma::Array{Int64,1})
	sigma_inv = zeros(Int64, length(sigma))
	for iii=1:length(sigma)
		sigma_inv[sigma[iii]] = iii
	end
	return sigma_inv
end

function num_cliques(sntree::SuperNodeTree)
	return length(sntree.snd_post)
end

# elimination tree algorithm from H.Liu - A Compact Row Storage Scheme for Cholesky Factors Using Elimination Trees
# function etree_liu(g)
# 	N = length(g.adjacency_list)
# 	par = zeros(Int64,N)
# 	ancestor = zeros(Int64,N)

# 	elemSequence = g.reverse_ordering[collect(1:N)]
# 	for iii in elemSequence
# 		for vk in g.adjacency_list[iii]
# 			if g.ordering[vk] < g.ordering[iii]
# 				r = vk
# 				while (ancestor[r] != 0) && (ancestor[r] != iii)
# 					t = ancestor[r]
# 					ancestor[r] = iii
# 					r = t
# 				end
# 				if ancestor[r] == 0
# 					ancestor[r] = iii
# 					par[r] = iii
# 				end
# 			end
# 		end
# 	end

# 	return par
# end

# simplified version of my own elimination tree algorithm with simplified data structure (fastest)
function etree(L)
	N = size(L, 1)
	par = zeros(Int64, N)
	# loop over Vertices of graph
	for i=1:N
		value = i
		# number of i-neighbors with order higher than order of node i
		par_ = find_parent_direct(L, i)
		par[i] = par_
	end
	return par
end

# perform a depth-first-search to determine the post order of the tree defined by parent and children vectors
# FIXME: This can be made a lot faster for the case that merges happened, i.e. Nc != length(par)

function post_order(par::Array{Int64,1}, child::Array{Array{Int64,1}}, Nc::Int64)

	order = (Nc + 1) * ones(Int64, length(par))
	root = findall(x ->x == 0, par)[1]
	stack = Array{Int64}(undef, 0)
	iii = Nc
	push!(stack, root)
	while !isempty(stack)
		v = pop!(stack)
		order[v] = iii
		iii -= 1
		push!(stack, child[v]...)
	end
	post = collect(1:length(par))
	sort!(post, by = x-> order[x])

	# if merges happened, remove the entries pointing to empty arrays / cliques
	Nc != length(par) && resize!(post, Nc)
	return post
end
post_order(par::Array{Int64,1}, child::Array{Array{Int64,1}}) = post_order(par, child, length(par))


function child_from_par(par::Array{Int64,1})

	child = [Array{Int64}(undef, 0) for i = 1:length(par)]
	for i = 1:length(par)
		par_ = par[i]
		par_ != 0 && push!(child[par_], i)
	end
	return child
end

# Using the post order ensures that no empty arrays from the clique merging are returned
function get_snd(sntree::SuperNodeTree, ind::Int64)
		return sntree.snd[sntree.snd_post[ind]]
end

function get_sep(sntree::SuperNodeTree, ind::Int64)
		return sntree.sep[sntree.snd_post[ind]]
end

# the block sizes are stored in post order, e.g. if clique 4 (stored in pos 4) has order 2, then nBlk[2] represents the size of clique 4
function get_nBlk(sntree::SuperNodeTree, ind::Int64)
		return sntree.nBlk[ind]
end

function get_clique(sntree::SuperNodeTree, ind::Int64)
	c = sntree.snd_post[ind]
	return union(sntree.snd[c], sntree.sep[c])
end

function print_cliques(sntree::SuperNodeTree)
	N = length(sntree.snd)

	println("Cliques of Graph:")
	for iii = 1:N
		println("$(iii): res: $(sntree.snd[iii]), sep: $(sntree.sep[iii])")
	end
end

function print_supernodes(sntree::SuperNodeTree)
	N = length(sntree.snd)
	println("Supernodes of Graph:")
	for iii = 1:N
		println("$(iii): $(sntree.snd[iii])")
	end
end


function check_degree_condition(v::Int64, w::Int64, degrees::Array{Int64,1})
	return degrees[v] > degrees[w] - 1
end


# Algorithm from A. Poten and C. Sun: Compact Clique Tree Data Structures in Sparse Matrix Factorizations (1989)
function pothen_sun(par::Array{Int64,1}, child::Array{Array{Int64,1}}, post::Array{Int64,1}, degrees::Array{Int64,1})
	N = length(par)
	sn_ind = -1 * ones(Int64, N) # if snInd[v] < 0 then v is a rep vertex, otherwise v ∈ supernode[snInd[v]]
	supernode_par = -1 * ones(Int64, N)

	# go through vertices in post_order
	for v in post
		child_ind = 0
		# check degree condition for all of v's childs from highest to lowest

		for (iii, w) in enumerate(child[v])
			# if not deg+(v) > deg+(w) - 1 for a certain w, set u to be w in snd(u), add v to snd(u)
			if !check_degree_condition(v, w, degrees)
				sn_ind[w] < 0 ? (u = w) : (u = sn_ind[w])
				sn_ind[v] = u
				child_ind = iii
				break
			end
		end

		# if v is still a rep vertex (i.e. above loop didnt find a child that fulfilled degree condition)
		if sn_ind[v] < 0
			u = v
		end

		for (iii, w) in enumerate(child[v])
			if iii != child_ind
				sn_ind[w] < 0 ? (x = w) : x = sn_ind[w]
				supernode_par[x] = u
			end
		end
	end

	# representative vertices
	reprv = findall(x-> x < 0, sn_ind)
	# vertices that are the parent of representative vertices
	repr_par = supernode_par[reprv]
	# take into account that all non-representative arrays are removed from the parent structure
	sn_par = zeros(Int64, length(reprv))

	for (iii, rp) in enumerate(repr_par)
		ind = findfirst(x -> x == rp, reprv)
		ind == nothing && (ind = 0)
		sn_par[iii] = ind
	end

	return sn_par, sn_ind
end

function find_supernodes(par::Array{Int64,1}, child::Array{Array{Int64,1}}, post::Array{Int64,1}, degrees::Array{Int64,1})
	supernode_par, snInd = pothen_sun(par, child, post, degrees)
	# number of vertices
	N = length(par)
	# number of representative vertices == number of supernodes
	Nrep = length(supernode_par)
	snode = [Array{Int64}(undef, 0) for i = 1:Nrep]

	for iii in post
		f = snInd[iii]
		if f < 0
			push!(snode[iii], iii)

		else
			push!(snode[f], iii)
		end
	end
	return snode, supernode_par

end

function find_separators(L, snodes::Array{Array{Int64,1},1}, supernode_par::Array{Int64,1}, post::Array{Int64,1})
	postInv = invert_order(post)

	Nc = length(supernode_par)
	sep = [Array{Int64}(undef, 0) for i = 1:Nc]

	for iii = 1:Nc
		vRep = snodes[iii][1]

		adjPlus = find_higher_order_neighbors(L, vRep)
		deg = length(adjPlus) + 1
		sep[iii] = adjPlus
		setdiff!(sep[iii], snodes[iii])

	end

	return sep

end

# Given a sparse lower triangular matrix L find the neighboring vertices of v
function find_neighbors(L::SparseMatrixCSC, v::Int64)
	L = Symmetric(L)
	col_ptr = L.colptr
	row_val = L.rowval
	return row_val[col_ptr[v]:col_ptr[v + 1] - 1]
end

function find_higher_order_neighbors(L::SparseMatrixCSC, v::Int64)
	v == size(L, 1) && return 0
	col_ptr = L.colptr
	row_val = L.rowval
	return row_val[col_ptr[v]:col_ptr[v + 1] - 1]
end


function find_parent_direct(L::SparseMatrixCSC, v::Int64)
	v == size(L, 1) && return 0
	col_ptr = L.colptr
	row_val = L.rowval
	return row_val[col_ptr[v]]
end



# findall the cardinality of adj+(v) for all v in V
function higher_degrees(L::SparseMatrixCSC)
	N = size(L, 1)
	degrees = zeros(Int64, N)
	col_ptr = L.colptr
	row_val = L.rowval

	for v = 1:N-1
		degrees[v] = col_ptr[v + 1] - col_ptr[v]
	end
	return degrees
end


# -------------------------------------
# FUNCTION DEFINITIONS
# -------------------------------------

# this assumes a sparse lower triangular matrix L
function connect_graph!(L::SparseMatrixCSC)
	# unconnected blocks don't have any entries below the diagonal in their right-most column
	m = size(L, 1)
	row_val = L.rowval
	col_ptr = L.colptr
	for j = 1:m-1
		connected = false
		for k in col_ptr[j]:col_ptr[j+1]-1
			if row_val[k] > j
				connected  = true
				break
			end
		end
		if !connected
			L[j+1, j] = 1
		end
	end
end

function find_graph!(ci, rows::Array{Int64, 1}, N::Int64, C::AbstractConvexSet)
	row_val, col_val = row_ind_to_matrix_indices(rows, N, C)
	F = QDLDL.qdldl(sparse(row_val, col_val, ones(length(row_val))), logical = true)
	# this takes care of the case that QDLDL returns an unconnected adjacency matrix L
	connect_graph!(F.L)
	ci.L = F.L
	return F.perm
end


# given an array [rows] that represent the nonzero entries of a vectorized NxN matrix,
# return the rows and columns of the nonzero entries of the original matrix
function row_ind_to_matrix_indices(rows::Array{Int64,1}, N::Int64, ::PsdCone{Float64})
	row_val = zeros(Int64, length(rows))
	col_val = zeros(Int64, length(rows))
	for (ind, r) in enumerate(rows)
		_rem = mod(r, N)
		fl = fld(r, N)
		if _rem == 0
			row_val[ind] = N
			col_val[ind] = fl
		else
			row_val[ind] = _rem
			col_val[ind] = fl + 1
		end
	end
	return row_val, col_val
end

# given an array [rows] that represent the nonzero entries of the vectorized upper triangular part of a NxN matrix,
# return the rows and columns of the nonzero entries of the original matrix
function row_ind_to_matrix_indices(rows::Array{Int64,1}, N::Int64, ::COSMO.PsdConeTriangle{Float64})
	#  allocate conservative here since we don't know how many diagonal entries are contained in row
	row_val = zeros(Int64, 2 * length(rows) + N)
	col_val = zeros(Int64, 2 * length(rows) + N)
	ind = 1
	_step = 1
	for (iii, r) in enumerate(rows)
		i, j = COSMO.svec_to_mat(r)

		row_val[ind] = i
		col_val[ind] = j
		ind += 1

		if i != j
			row_val[ind] = j
			col_val[ind] = i
			ind += 1
		end
	end
	return row_val[1:ind - 1], col_val[1:ind - 1]
end

# Given a linear index find the col of the corresponding upper triangular matrix
function svec_get_col(x::Int64)
	c = (sqrt(8 * x + 1) - 1) / 2
	if c % 1 != 0
		return Int((floor(c) + 1))
	else
		return Int(c)
	end
end

# Given a linear index find the row of the corresponding upper triangular matrix
function svec_get_row(x::Int64)
	c = get_col(x) - 1
	k = (c + 1) * c / 2
	return (x - k)::Int64
end

# Given a linear index find the row and column of the corresponding upper triangular matrix
function svec_to_mat(ind::Int64)
	c = svec_get_col(ind)
	k = (c - 1) * c / 2
	r = Int(ind - k)
	return r::Int64, c::Int64
end

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
    NoMerge() <: AbstractMergeStrategy
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
		@show(cand)

		# evaluate wether to merge the candidates
		if evaluate(t, strategy, cand)
			merge_cliques!(t, cand)
		end

		strategy.stop && break
		# update strategy information after the merge
		update!(strategy, t, cand)
	end

	# number of cliques
	Nc = length(findall(x -> x != -1, t.snd_par))
	# the merging creates empty supernodes and seperators, recalculate a post order for the supernodes
 	 snd_post = post_order(t.snd_par, t.snd_child, Nc)
 	 t.snd_post = snd_post

 	# recalculate block sizes (attention: the block sizes are stored in post order)
	t.nBlk = zeros(Nc)
  for iii = 1:Nc
  	c = snd_post[iii]
		t.nBlk[iii] = Base.power_by_squaring(length(t.sep[c]) + length(t.snd[c]), 2)
  end
  return nothing
end


function merge_cliques!(t::SuperNodeTree, cand::Tuple{Int64,Int64})
	cand[1] < cand[2] ? merge_child!(t, cand) : merge_sibling!(t, cand)
end

# Merge two cliques that are in a parent (cand[1]) - child (cand[2]) relationship
function merge_child!(t::SuperNodeTree, cand::Tuple{Int64,Int64})
	p = cand[1]
	ch = cand[2]

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

	return nothing
end

# Merge two cliques that are in a sibling relationship
function merge_sibling!(t::SuperNodeTree, cand::Tuple{Int64,Int64})
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

	return nothing
end

"""
    NoMerge() <: AbstractMergeStrategy
A strategy that does not merge cliques.
"""
 struct NoMerge <: AbstractMergeStrategy end


"""
    PairwiseMerge <: AbstractMergeStrategy

A merge strategy that calculates the edge metric `A ∩ B / A ∪ B` for every two cliques that are in a parent-child or sibling relationship. The resulting clique
graph is traversed from the highest edge metric to the lowest.
"""
mutable struct PairwiseMerge <: AbstractMergeStrategy
	stop::Bool
	edges::AbstractMatrix
	edge_score::AbstractEdgeScore

	function PairwiseMerge(; edge_score = RelIntersect)
		new(false, spzeros(Float64, 0, 0), edge_score)
	end
end

# compute the edge set of the initial clique graph, only consider parent-child and sibling-sibling relationships
function initialise!(t, strategy::PairwiseMerge)
	n_cliques = length(t.snd_par)
	strategy.edges = spzeros(Float64, n_cliques, n_cliques)

	# following a descending post-order, add edges
	for i = 1:length(t.snd_post)
		clique_ind = t.snd_post[i]
		children = t.snd_child[clique_ind]
		n_children = length(children)
		# add edge to each child (father: row, child: col) (Upper Triangle)
		for (j, ch) in enumerate(children)
			strategy.edges[clique_ind, ch] = edge_metric_parent(t.snd, t.sep, clique_ind, ch, strategy.edge_score)

			# add edges between siblings combinations (Lower Triangle)
			for sib in view(children,(j+1):n_children)
				strategy.edges[sib, ch] = edge_metric_siblings(t.snd, t.sep, ch, sib, strategy.edge_score)
			end
		end
	end
	return nothing
end

traverse(t, strategy::PairwiseMerge) = traverse(t, strategy, strategy.edge_score)
function traverse(t, strategy::PairwiseMerge, edge_score::RelIntersect)
	~, ind = findmax(strategy.edges)
	return (ind[1], ind[2])
end

function traverse(t, strategy::PairwiseMerge, edge_score::ComplexityScore)
	~, ind = findmax(strategy.edges.nzval)
	r = strategy.edges.rowval[ind]
	c = findfirst(x -> x >= ind, strategy.edges.colptr)
	if strategy.edges.colptr[c] > ind
		c -= 1
	end
	return (r, c)
end

evaluate(t, strategy::PairwiseMerge, cand) = evaluate(t, strategy, cand, strategy.edge_score)

function evaluate(t, strategy::PairwiseMerge, cand, edge_score::RelIntersect)
	c1 = cand[1]
	c2 = cand[2]

	# individual block sizes
	n_1 = length(t.snd[c1]) + length(t.sep[c1])
	n_2 = length(t.snd[c2]) + length(t.sep[c2])

	# merged block size
	n_m = n_1 + length(t.snd[c2])
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
compute_complexity_savings(n_1::Int64, n_2::Int64, n_m::Int64) = n_1^3 - n_2^3 - n_m^3

# Approximates the number of operations for one projection of all the cliques in the tree
compute_complexity(t::COSMO.SuperNodeTree) = sum(map(x -> x^3, t.nBlk))



# After a merge operation update the information of the strategy
function update!(strategy::PairwiseMerge, t, cand)
	c1 = cand[1]
	c2 = cand[2]
	edges = strategy.edges

	# remove edge value between c1 and c2 and all to c2
	edges[c1, c2] = 0
	@. edges[c2, :] = 0
	@. edges[:, c2] = 0

	# recalculate edge values of c1's children
	for (j, ch) in enumerate(t.snd_child[c1])
		edges[c1, ch] = COSMO.edge_metric_parent(t.snd, t.sep, c1, ch)
	end

	# edge to c1's parent (for all non-root nodes)
	p = t.snd_par[c1]
	if p > 0
		edges[p, c1] = COSMO.edge_metric_parent(t.snd, t.sep, p, c1)

		# edges to c1's siblings
		siblings = t.snd_child[p]
		for sib in siblings
			if sib != c1
				# make sure it's stored in the lower triangle of matrix edges
				if sib < c1
					edges[c1, sib] = COSMO.edge_metric_siblings(t.snd, t.sep, c1, sib)
				else
					edges[sib, c1] = COSMO.edge_metric_siblings(t.snd, t.sep, c1, sib)
				end
			end
		end
	end
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
	# individual block sizes
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
	n_m = n_1 + n_2 - intersect_dim(sep[c_a], sep[c_b])
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

initialise!(t, strategy::TreeTraversalMerge) = nothing

# traverse tree in topological order and return clique and its parent
function traverse(t, strategy::TreeTraversalMerge)
	c = t.snd_post[strategy.clique_ind]
	return (t.snd_par[c], c)
end

function evaluate(t, strategy::TreeTraversalMerge, cand)
	c = cand[1]
	par = cand[2]
	n_sep = length(t.sep[c])
	n_sep_par = length(t.sep[par])
	n_cl = length(t.snd[c]) + n_sep
	n_cl_par = length(t.snd[par]) + n_sep
	return (n_cl_par - n_sep) * (n_cl - n_sep) <= strategy.t_fill || max(n_cl - n_sep, n_cl_par - n_sep_par) <= strategy.t_size
end

function update!(strategy::TreeTraversalMerge, t, cand)
	# break if node was last one
	if strategy.clique_ind == length(t.snd)
		strategy.stop = true
	# otherwise increment node index
	else
		strategy.clique_ind += 1
	end
end