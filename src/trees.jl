export print_merge_logs, print_clique_sizes

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



abstract type AbstractTreeBasedMerge <: AbstractMergeStrategy end
abstract type AbstractGraphBasedMerge <: AbstractMergeStrategy end


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
	strategy::AbstractMergeStrategy
	function SuperNodeTree(L::SparseMatrixCSC,  merge_strategy::AbstractMergeStrategy)
		par = etree(L)
		child = child_from_par(par)
		post = post_order(par, child)
		# sorting of children from highest order one to lowest make sure that a node always gets
		# added to the supernode with the highest rep vertex
		#sort_children!(child, post)
		# faster algorithm according to Vandenberghe p.308
		degrees = higher_degrees(L)
		snd, snd_par = find_supernodes(par, post, degrees)

		snd_child = child_from_par(snd_par)
	 	snd_post = post_order(snd_par, snd_child)

		# If the merge strategy is clique graph-based, we give up the tree structure and add the seperators to the supernodes
		# the supernodes then represent the full clique
		# after the clique merging a new clique tree will be computed before psd completion is performed
	 	if typeof(merge_strategy) <: AbstractGraphBasedMerge
			add_separators!(L, snd, snd_par, post)
			@. snd_par = -1
			snd_child = [Int64[] for i = 1:length(snd)]
			sep = [Int64[] for i = 1:length(snd)]
			new(snd, snd_par, snd_post, snd_child, post, par, sep, [0], length(snd_post), MergeLog(), merge_strategy)
	 	# If the merge strategy is tree based keep the supernodes and separators in to different locations
		else
	 		# given the supernodes (clique residuals) find the clique separators
			sep = find_separators(L, snd, snd_par, post)
			new(snd, snd_par, snd_post, snd_child, post, par, sep, [0], length(snd_post), MergeLog(), merge_strategy)

		end
	end


	# FIXME: only for debugging purposes
	function SuperNodeTree(res, par, snd_post, sep, merge_strategy; post::Array{Int64, 1} = [1])
		child = child_from_par(par)
  	new(res, par, snd_post, child, post, [1], sep, [1], length(res), MergeLog(), merge_strategy)
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

# ------------------------------------------
# Interfaces function to SuperNodeTree
# ------------------------------------------

function num_cliques(sntree::SuperNodeTree)
	return sntree.num::Int64
end

function get_post_order(sntree::SuperNodeTree)
	return sntree.snd_post
end

function get_post_order(sntree::SuperNodeTree, i::Int64)
	return sntree.snd_post[i]
end
# Using the post order ensures that no empty arrays from the clique merging are returned
function get_snd(sntree::SuperNodeTree, ind::Int64)
		return sntree.snd[sntree.snd_post[ind]]
end

function get_sep(sntree::SuperNodeTree, ind::Int64)
		return sntree.sep[sntree.snd_post[ind]]
end

# the block sizes are stored in post order, e.g. if clique 4 (stored in pos 4) has order 2, then nBlk[2] represents the cardinality of clique 4
function get_nBlk(sntree::SuperNodeTree, ind::Int64)
		return sntree.nBlk[ind]::Int64
end

"Returns the number of rows of all the blocks (cliques) represented in the tree after decomposition."
function get_decomposed_dim(sntree::SuperNodeTree, C::DecomposableCones{<: Real})
	dim = 0
	for iii = 1:num_cliques(sntree)
		dim += vec_dim(get_nBlk(sntree, iii), C)
	end
	return dim::Int64
end

"Given the side dimension of a PSD cone return the number of stored entries."
vec_dim(side_dim::Int64, C::PsdCone{<:Real}) = Base.power_by_squaring(side_dim, 2)
vec_dim(side_dim::Int64, C::PsdConeTriangle{<:Real}) = div(side_dim * (side_dim + 1), 2)

" Return clique with post order `ind` (prevents returning empty arrays due to clique merging)"
get_clique(sntree::SuperNodeTree, ind::Int64) = get_clique(sntree, ind, sntree.strategy)
function get_clique(sntree::SuperNodeTree, ind::Int64, strategy::AbstractTreeBasedMerge)
	c = sntree.snd_post[ind]
	return union(sntree.snd[c], sntree.sep[c])
end

function get_clique(sntree::SuperNodeTree, ind::Int64, strategy::AbstractGraphBasedMerge)
	# if a clique tree has been recomputed, call the method for AbstractTreeBasedMerge types
	if strategy.clique_tree_recomputed
		return get_clique(sntree, ind, COSMO.NoMerge())
	else
		c = sntree.snd_post[ind]
		return sntree.snd[c]
	end
end


function print_cliques(sp; reordered = true)
	sntree = sp.sntree
	reverse_ordering = sp.reverse_ordering
	Nsnd = length(sntree.snd)
	println("Cliques of Graph:")
	println("Reordered = $(reordered)")
	for iii = 1:Nsnd
			if !reordered
				snd = map(x-> reverse_ordering[x], sntree.snd[iii])
				sep = map(x-> reverse_ordering[x], sntree.sep[iii])
			else
				snd = sntree.snd[iii]
				sep = sntree.sep[iii]
			end
			println("$(iii): res: $(snd), sep: $(sep)")
	end
end
print_cliques(sntree::SuperNodeTree) = print_cliques(sntree, sntree.strategy)
function print_cliques(sntree::SuperNodeTree, strategy::AbstractTreeBasedMerge)
	Nsnd = length(sntree.snd)
	println("Cliques of Graph:")
	for iii = 1:Nsnd
			println("$(iii): res: $(sntree.snd[iii]), sep: $(sntree.sep[iii])")
	end
end

function print_clique_orders(sntree::SuperNodeTree, strategy::AbstractTreeBasedMerge)
	Nsnd = length(sntree.snd)
	println("Cliques of Graph:")
	for iii = 1:Nsnd
		post = sntree.post
		postInv = invert_order(post) #σ-1(v) = j

  	snd_ip = map(x -> postInv[x], sntree.snd[iii])
  	sep_ip = map(x -> postInv[x], sntree.sep[iii])
		println("$(iii): res: $(snd_ip), sep: $(sep_ip)")
	end
end

function print_cliques(sntree::SuperNodeTree, strategy::AbstractGraphBasedMerge)

	if strategy.clique_tree_recomputed
		return print_cliques(sntree, COSMO.NoMerge())
	else
		Nsnd = length(sntree.snd)
		println("Cliques of Graph:")
		for iii = 1:Nsnd
				println("$(iii): clique: $(sntree.snd[iii])")
		end
	end
	return nothing
end


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



function check_degree_condition(v::Int64, w::Int64, degrees::Array{Int64,1})
	return degrees[v] > degrees[w] - 1
end


# Algorithm from A. Poten and C. Sun: Compact Clique Tree Data Structures in Sparse Matrix Factorizations (1989)
function pothen_sun(par::Array{Int64,1}, post::Array{Int64,1}, degrees::Array{Int64,1})
	N = length(par)
	sn_ind = -1 * ones(Int64, N) # if snInd[v] < 0 then v is a rep vertex, otherwise v ∈ supernode[snInd[v]]
	supernode_par = -1 * ones(Int64, N)
	children = 	[Array{Int64}(undef, 0) for i = 1:length(par)]

	root_ind = findfirst(x -> x == 0, par)
	# go through parents of vertices in post_order
	for v in post

		if par[v] == 0
			push!(children[root_ind], v)
		else
			push!(children[par[v]], v)
		end

		# parent is not the root
		if par[v] != 0
			if degrees[v] - 1 == degrees[par[v]] && sn_ind[par[v]] == -1
				# Case A: v is a representative vertex
				if sn_ind[v] < 0
					sn_ind[par[v]] = v
					sn_ind[v] -= 1
				# Case B: v is not representative vertex, add to sn_ind[v] instead
				else
					sn_ind[par[v]] = sn_ind[v]
					sn_ind[sn_ind[v]] -= 1
				end
			else
				if sn_ind[v] < 0
					supernode_par[v] = v
				else
					supernode_par[sn_ind[v]] = sn_ind[v]
				end
			end
		end

		# k: rep vertex of the snd that v belongs to
		if sn_ind[v] < 0
			k = v
		else
			k = sn_ind[v]
		end
		# loop over v's children
		v_children = children[v]
		if !isempty(v_children)
			for w in v_children
				if sn_ind[w] < 0
					l = w
				else
					l = sn_ind[w]
				end
				if l != k
					supernode_par[l] = k
				end
			end
		end
	end # loop over vertices

	# 	for (iii, w) in enumerate(child[v])
	# 		for
	# 	end
	# 	child_ind = 0
	# 	# check degree condition for all of v's children from highest to lowest

	# 	for (iii, w) in enumerate(child[v])
	# 		# if not deg+(v) > deg+(w) - 1 for a certain w, set u to be w in snd(u), add v to snd(u)
	# 		if !check_degree_condition(v, w, degrees)
	# 			sn_ind[w] < 0 ? (u = w) : (u = sn_ind[w])
	# 			sn_ind[v] = u
	# 			child_ind = iii
	# 			break
	# 		end
	# 	end

	# 	# if v is still a rep vertex (i.e. above loop didnt find a child that fulfilled degree condition)
	# 	if sn_ind[v] < 0
	# 		u = v
	# 	end

	# 	for (iii, w) in enumerate(child[v])
	# 		if iii != child_ind
	# 			sn_ind[w] < 0 ? (x = w) : x = sn_ind[w]
	# 			supernode_par[x] = u
	# 		end
	# 	end
	# end

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

function find_supernodes(par::Array{Int64,1}, post::Array{Int64,1}, degrees::Array{Int64,1})
	supernode_par, snInd = pothen_sun(par, post, degrees)
	# number of vertices
	N = length(par)
	# number of representative vertices == number of supernodes
	Nrep = length(supernode_par)
	snode = [Array{Int64}(undef, 0) for i = 1:N]

	for iii = 1:N
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

function add_separators!(L::SparseMatrixCSC, snodes::Array{Array{Int64,1},1}, supernode_par::Array{Int64,1}, post::Array{Int64,1})
	postInv = invert_order(post)

	Nc = length(supernode_par)

	for iii = 1:Nc
		snode = snodes[iii]
		vRep = snode[1]

		adjPlus = find_higher_order_neighbors(L, vRep)
		for neighbor in adjPlus
			if !in(neighbor, snode)
				push!(snode, neighbor)
			end
		end
		sort!(snode)
	end
	return nothing
end


"""
		reorder_snd_consecutively!(sntree, ordering)

Takes a SuperNodeTree and reorders the vertices in each supernode (and separator) to have consecutive order.

The reordering is needed to achieve equal column structure for the psd completion of the dual variable `Y`. This also modifies `ordering` which maps the vertices in the sntree back to the actual location in the not reordered data, i.e.
the primal constraint variable `S` and dual variables `Y`.
"""
function reorder_snd_consecutively!(t::SuperNodeTree, ordering::Array{Int64, 1})

	# determine permutation vector p and permute the vertices in each snd
	p = zeros(Int64,length(t.post))
	snd = t.snd
	sep = t.sep

	k = 1
	for snd_ind in t.snd_post
	  s = snd[snd_ind]
	  l = length(s)
	  p[k:k+l-1] = s
	  s[:] = collect(k:1:k+l-1)
	  k += l
	end


  # permute the separators as well
  p_inv = invperm(p)
  for sp in t.sep
    sp[:] = map(x -> p_inv[x], sp)
  end

	# # permute the degrees, supernodal parent structure and supernodal post order

	# snd_post_inv = invperm(t.snd_post)
	# root_ind = findfirst(x -> x == 0, t.snd_par)
	# t.snd_par[root_ind] = root_ind
	# t.snd_par = sndpost_inv[t.snd_par[t.snd_post]]
	# t.snd_post = collect(1:1:t.num)

	permute!(ordering, p)
	return nothing
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
