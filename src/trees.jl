export print_merge_logs, print_clique_sizes

"""
		MergeLog

A struct to analyse the clique merges. Introduced for debugging purposes.
"""
mutable struct MergeLog
	num::Int64 # number of merges
	clique_pairs::Array{Int64, 2} # ordered pair merges
	decisions::Array{Bool, 1} # at what step was merged
	function MergeLog()
		new(0, zeros(0, 2), Array{Bool}(undef, 0))
	end
end

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

function print_clique_sizes(ws)
	sp_arr = ws.ci.sp_arr
	# print the merge logs for each sparsity pattern
	println(">>> Clique Dimensions:")
	for (iii, sp) in enumerate(sp_arr)
		println("Sparsity Pattern Nr. $(iii), Graph Size: $(length(sp.sntree.par))")

		t = sp.sntree
		Nc = length(t.snd_post)
		# block sizes
		sizes = zeros(Int64, Nc)
		# occurences of block size
		occ = zeros(Int64, Nc)
		for jjj = 1:Nc
			c = t.snd_post[jjj]
			dim = length(t.snd[c]) + length(t.sep[c])
			# try to check if that dimension occured before, if not add, otherwise increment occ
			ind = findfirst(x -> x == dim, sizes)

			if ind != nothing
				occ[ind] = occ[ind] + 1
			else
				sizes[jjj] = dim
				occ[jjj] = 1
			end
		end

		# consolidate
		filter!(x -> x != 0, occ)
		filter!(x -> x != 0, sizes)

		# sort sizes
		p = sortperm(sizes)
		sizes = sizes[p]
		occ = occ[p]

		for (jjj, dim) in enumerate(sizes)
			println("$(occ[jjj])x dim:$(dim)")
		end
	end
end
"""
		SuperNodeTree

A structure to represent and analyse the sparsity pattern of the input matrix L.
"""
mutable struct SuperNodeTree
	snd::Array{Array{Int64,1},1} #vertices of supernodes stored in one array (also called residuals)
	snd_par::Array{Int64,1}  # parent of supernode k is supernode j=snd_par[k]
	snd_post::Array{Int64,1} # post order of supernodal elimination tree
	snd_child::Array{Array{Int64,1},1}
	post::Array{Int64} # post ordering of the vertices in elim tree σ(j) = v
	par::Array{Int64}
	sep::Array{Array{Int64,1},1} #vertices of clique seperators
	nBlk::Array{Int64,1} #sizes of submatrizes defined by each clique, sorted by post-ordering, e.g. size of clique with order 3 => nBlk[3]
	num::Int64 # number of supernodes / cliques in tree
	merge_log::MergeLog
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

		new(snd, snd_par, snd_post, snd_child, post, par, sep, [0], length(snd_post), MergeLog())
	end
	# FIXME: only for debugging purposes
	function SuperNodeTree(res, par, post, sep)
		child = child_from_par(par)
  	new(res, par, post, child, [1], [1], sep, [1], length(res), MergeLog())
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
	return sntree.num
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
	snode = [Array{Int64}(undef, 0) for i = 1:N]

	for iii in post
		f = snInd[iii]
		if f < 0
			push!(snode[iii], iii)

		else
			push!(snode[f], iii)
		end
	end
	filter!(x -> !isempty(x), snode)
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
	row_val, col_val = COSMO.row_ind_to_matrix_indices(rows, N, C)
	F = QDLDL.qdldl(sparse(row_val, col_val, ones(length(row_val))), logical = true)#, perm = collect(1:N))
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
    NoMerge <: AbstractMergeStrategy

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


function update!(strategy::TreeTraversalMerge, t, cand)
	# try to merge last node of order 1, then stop
	if strategy.clique_ind == 1
		strategy.stop = true
	# otherwise decrement node index
	else
		strategy.clique_ind -= 1
	end
end