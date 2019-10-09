function _contains(convex_sets::COSMO.CompositeConvexSet, t::Type{<:COSMO.AbstractConvexCone})
  for set in convex_sets.sets
    if typeof(set) <: t
      return true
    end
  end
  return false
end

function chordal_decomposition!(ws::COSMO.Workspace)
  # do nothing if no psd cones present in the problem
  if !_contains(ws.p.C, DecomposableCones{Float64})
    ws.settings.decompose = false
    return nothing
  end
  ws.ci = ChordalInfo{Float64}(ws.p)

  find_sparsity_patterns!(ws)

  if ws.ci.num_decomposable > 0

    # Do transformation similar to clique tree based transformation in SparseCoLo
    if ws.settings.colo_transformation

      augment_clique_based!(ws)



    else
      # find transformation matrix H and new composite convex set
      find_decomposition_matrix!(ws)

      # augment the original system
      augment_system!(ws)
    end
    pre_allocate_variables!(ws)

  else
    ws.settings.decompose = false
  end
  nothing
end

# analyse PSD cone constraints for chordal sparsity pattern
function find_sparsity_patterns!(ws)
  row_ranges = get_set_indices(ws.p.C.sets)
  sp_ind = 1
  for (k, C) in enumerate(ws.p.C.sets)
    psd_row_range = row_ranges[k]
    csp = find_aggregate_sparsity(ws.p.A, ws.p.b, psd_row_range, C)
    sp_ind = analyse_sparsity_pattern!(ws.ci, csp, ws.p.C.sets, C, k, psd_row_range, sp_ind, ws.settings.merge_strategy)
  end
end

analyse_sparsity_pattern!(ci, csp, sets, C::AbstractConvexSet, k, psd_row_range, sp_ind, merge_strategy) = sp_ind

function analyse_sparsity_pattern!(ci, csp, sets, C::DecomposableCones{T}, k, psd_row_range, sp_ind, merge_strategy) where {T <: Real}
  if length(csp) < C.dim
    return _analyse_sparsity_pattern(ci, csp, sets, C, k, psd_row_range, sp_ind, merge_strategy)
  else
   sets[k] = COSMO.DenseEquivalent(C, C.dim)
   return sp_ind
 end
end

function _analyse_sparsity_pattern(ci, csp, sets, C::Union{PsdCone{<: Real}, PsdConeTriangle{<: Real}}, k, psd_row_range, sp_ind, merge_strategy) where {T <: Real}
  ordering, nz_ind_map = find_graph!(ci, csp, C.sqrt_dim, C)
  sp = COSMO.SparsityPattern(ci.L, C.sqrt_dim, ordering, merge_strategy, psd_row_range, k, nz_ind_map)
  # if after analysis of SparsityPattern & clique merging only one clique remains, don't bother decomposing
  if num_cliques(sp.sntree) == 1
    sets[k] = DenseEquivalent(C, C.dim)
    return sp_ind
  else
    ci.sp_arr[sp_ind] = sp
    push!(ci.psd_cones_ind, k)
    ci.num_decomposable += 1
    return sp_ind + 1
  end
end

DenseEquivalent(C::COSMO.PsdCone{T}, dim::Int64) where {T} = COSMO.DensePsdCone{T}(dim)
DenseEquivalent(C::COSMO.PsdConeTriangle{T}, dim::Int64) where {T} = COSMO.DensePsdConeTriangle{T}(dim)




function nz_rows(a::SparseMatrixCSC, ind::UnitRange{Int64}, DROP_ZEROS_FLAG::Bool)
  DROP_ZEROS_FLAG && dropzeros!(a)
  active = falses(length(ind))
  for r in a.rowval
    if in(r, ind)
      active[r - ind.start + 1] = true
    end
  end
  return findall(active)
end

function number_of_overlaps_in_rows(A::SparseMatrixCSC)
  # sum the entries row-wise
  numOverlaps = sum(A, dims = 2)
  ri = findall(x -> x > 1, numOverlaps)
  return ri, numOverlaps[ri]
end


function find_aggregate_sparsity(A, b, ind, C::DecomposableCones{ <: Real})
  AInd = nz_rows(A, ind, false)
  # commonZeros = AInd[find(x->x==0,b[AInd])]
  bInd = findall(x -> x != 0, view(b, ind))
  commonNZeros = union(AInd, bInd)
  return commonNZeros
end
find_aggregate_sparsity(A, b, ind, C::AbstractConvexSet) = Int64[]

function vec_to_mat_ind(ind::Int64, n::Int64)
  ind > n^2 && error("Index ind out of range.")
  ind == 1 && (return 1, 1)

  r = ind % n

  if r == 0
    j = Int(ind / n)
    i = n
  else
    j = Int(floor(ind / n) + 1)
    i = r
  end
  return i, j
end

# Converts the matrix element (i, j) of A ∈ m x n into the corresponding linear index of v = vec(A)
function mat_to_vec_ind(i::Int64, j::Int64, m::Int64)
  (i > m || i <= 0 || j <= 0) && throw(BoundsError("Indices outside matrix bounds."))
  return (j - 1) * m + i
end

# Converts the matrix element (i, j) of A ∈ m x n into the corresponding linear index of v = svec(A, ::UpperTriangle)
function mat_to_svec_ind(i::Int64, j::Int64)
  if i <= j
    return div((j - 1) * j, 2) + i
  else
    return div((i - 1) * i, 2) + j
  end
end

vectorized_ind(i::Int64, j::Int64, m::Int64, C::PsdCone{T}) where {T} = mat_to_vec_ind(i, j, m)
vectorized_ind(i::Int64, j::Int64, m::Int64, C::PsdConeTriangle{T}) where {T} = mat_to_svec_ind(i, j)

function svec_to_mat_ind(k::Int64)
  j = isqrt(2 * k)
  i = k - div((j - 1) * j, 2)
  return i, j
end


# function finds the transformation matrix H to decompose the vector s into its parts and stacks them into sbar
function find_decomposition_matrix!(ws)

  # allocate H and new decomposed cones
  n = COSMO.find_H_col_dimension(ws.p.C.sets, ws.ci.sp_arr)
  H_I = zeros(Int64, n)


  # find number of decomposed and total sets and allocate structure for new compositve convex set
  num_total, num_new_psd_cones = COSMO.num_cone_decomposition(ws)
  # decomposed_psd_cones = Array{COSMO.PsdCone}(undef, 0)
  C_new = Array{COSMO.AbstractConvexSet{Float64}}(undef, num_total)
  C_new[1] = COSMO.ZeroSet{Float64}(ws.ci.originalM)

  # loop over all convex sets and fill H and composite convex set accordingly
  row = 1
  entry = 1
  sp_ind = 1
  set_ind = 2
  for (kkk, C) in enumerate(ws.p.C.sets)
    set_ind, sp_ind, entry = COSMO.decompose!(H_I, C_new, set_ind, C, entry, row, ws.ci.sp_arr, sp_ind)
    row += C.dim
  end
  ws.p.C = COSMO.CompositeConvexSet(C_new)
  ws.ci.H = sparse(H_I, collect(1:n), ones(n))
end


function decompose!(H_I::Vector{Int64}, C_new, set_ind::Int64, C::COSMO.AbstractConvexSet, entry::Int64, row::Int64, sp_arr::Array{SparsityPattern}, sp_ind::Int64)
  #H[row:row + C.dim - 1, col:col + C.dim - 1] = sparse(1.0I, C.dim, C.dim)
  for i = 1:C.dim
    H_I[entry] = row + i - 1
    entry += 1
  end
  C_new[set_ind] = C

  return set_ind + 1, sp_ind, entry
end

# for a clique with nBlk entries, return the number of entries in the corresponding matrix
get_blk_rows(nBlk::Int64, C::COSMO.PsdCone{T}) where {T} = Base.power_by_squaring(nBlk, 2)
get_blk_rows(nBlk::Int64, C::COSMO.PsdConeTriangle{T}) where {T} = div((nBlk + 1) * nBlk, 2)

function decompose!(H_I::Vector{Int64}, C_new, set_ind::Int64,  C::DecomposableCones{ <: Real}, entry::Int64, row_start::Int64, sp_arr::Array{SparsityPattern}, sp_ind::Int64)
  sparsity_pattern = sp_arr[sp_ind]
  sntree = sparsity_pattern.sntree
  # original matrix size
  original_size = C.sqrt_dim

  for iii = 1:num_cliques(sntree)

    c = COSMO.get_clique(sntree, iii)
    # the graph and tree algorithms determined the cliques vertices of an AMD-permuted matrix. Since the location of the data hasn't changed in reality, we have to map the clique vertices back
    c = map(v -> sparsity_pattern.ordering[v], c)
    sort!(c)
    entry = COSMO.add_subblock_map!(H_I, c, original_size, entry, row_start, C)

    # create and add new cone for subblock
    num_rows = get_blk_rows(get_nBlk(sntree, iii), C)
    C_new[set_ind] = typeof(C)(num_rows)
    set_ind += 1
  end
  return set_ind, sp_ind + 1, entry
end



# fills the corresponding entries of H for clique c
function add_subblock_map!(H_I::Vector{Int64}, clique_vertices::Array{Int64}, original_size::Int64, entry::Int64, row_start::Int64, ::PsdCone{<: Real})
  for vj in clique_vertices
    for vi in clique_vertices
      row = mat_to_vec_ind(vi, vj, original_size)
      H_I[entry] = row_start + row - 1
      entry += 1
    end
  end
  return entry::Int64
end

function add_subblock_map!(H_I::Vector{Int64}, clique_vertices::Array{Int64}, original_size::Int64, entry::Int64,  row_start::Int64, ::PsdConeTriangle{<: Real})
  for vj in clique_vertices
    for vi in clique_vertices
      if vi <= vj
        row = mat_to_svec_ind(vi, vj)
        H_I[entry] = row_start + row - 1
        entry += 1
      end
    end
  end
  return entry::Int64
end


function find_H_col_dimension(sets, sp_arr)
  sum_cols = 0
  sp_arr_ind = 1
  for C in sets
    dim, sp_arr_ind = decomposed_dim(C, sp_arr, sp_arr_ind)
    sum_cols += dim
  end
  return sum_cols::Int64
end

"Return the dimension of the problem after a clique tree based decomposition, given the sparsity patterns in `sp_arr`."
function find_A_dimension(n_original::Int64, sets, sp_arr)
  num_cols = n_original
  num_overlapping_entries = 0
  num_rows = 0
  sp_arr_ind = 1
  for C in sets
    dim, overlaps, sp_arr_ind = decomposed_dim(C, sp_arr, sp_arr_ind)
    num_rows += dim
    num_overlapping_entries += overlaps
  end
  return num_rows::Int64, (num_cols + num_overlapping_entries)::Int64, num_overlapping_entries
end

decomposed_dim(C::AbstractConvexSet, sp_arr::Array{SparsityPattern}, sp_arr_ind::Int64) = (C.dim, 0, sp_arr_ind)
function decomposed_dim(C::DecomposableCones{ <: Real}, sp_arr::Array{SparsityPattern}, sp_arr_ind::Int64)
  sntree = sp_arr[sp_arr_ind].sntree
  dim, overlaps = get_decomposed_dim(sntree, C)
  return dim::Int64, overlaps::Int64, (sp_arr_ind + 1)::Int64
end


function num_cone_decomposition(ws)
  num_sets = length(ws.p.C.sets)
  num_old_psd_cones = length(ws.ci.psd_cones_ind)
  num_new_psd_cones = 0
  for iii = 1:num_old_psd_cones
    sp = ws.ci.sp_arr[iii]
    num_new_psd_cones += COSMO.num_cliques(sp.sntree)
  end
  ws.ci.num_decom_psd_cones = ws.ci.num_psd_cones - ws.ci.num_decomposable + num_new_psd_cones
  # the total number is the number of original non-psd cones + number of decomposed psd cones + 1 zeroset for the augmented system
  num_total = num_sets - num_old_psd_cones + num_new_psd_cones + 1
  return num_total, num_new_psd_cones
end

function augment_system!(ws)
  _augment!(ws.p, ws.ci.H)
  nothing
end

function _augment!(problem, H::SparseMatrixCSC)
  mH, nH = size(H)
  m, n = size(problem.A)
  problem.P = blockdiag(problem.P, spzeros(nH, nH))
  problem.q = vec([problem.q; zeros(nH)])
  problem.A = [problem.A H; spzeros(nH, n) -sparse(1.0I, nH, nH)]
  problem.b = vec([problem.b; zeros(nH)])
  problem.model_size[1] = size(problem.A, 1)
  problem.model_size[2] = size(problem.A, 2)
  nothing
end

function get_rows(b::SparseVector, row_range::UnitRange{Int64})
  rows = b.nzind
  if length(rows) > 0
    s = searchsortedfirst(rows, row_range.start)
    if rows[s] > row_range.stop || s == 0
        return nothing, 0:0
    else
      e = searchsortedlast(rows, row_range.stop)
      return s:e
    end
  else
    return nothing
  end

end

function get_rows(A::SparseMatrixCSC, col::Int64, row_range::UnitRange{Int64})
  erange = A.colptr[col]:(A.colptr[col + 1]-1)

  # if the column has entries
  if erange.start <= erange.stop
    # create a view into the row values of column col
    rows = view(A.rowval, erange)
    # find the rows within row_start:row_start+C.dim-1
    # s: index of first entry in rows >= row_start
    s = searchsortedfirst(rows, row_range.start)
    if rows[s] > row_range.stop || s == 0
      return 0:0
    else
      # e: index of last value in rows <= row_start + C.dim - 1
      e = searchsortedlast(rows, row_range.stop)
      return erange[s]:erange[e]
    end
  else
    return 0:0
  end
end


function add_entries!(A_I::Array{Int64, 1}, b_I::Array{Int64, 1}, C_new::Array{COSMO.AbstractConvexSet{Float64}, 1}, row_ptr::Int64, A0::SparseMatrixCSC, b0::SparseVector, row_range::UnitRange{Int64}, overlap_ptr::Int64, set_ind::Int64, sp_ind::Int64,
  sp_arr::Array{SparsityPattern, 1}, C::AbstractConvexSet, k::Int64, cone_map::Dict{Int64, Int64})

  m, n = size(A0)

  # let's handle vector b first
  offset = row_ptr - row_range.start

  row_range_col = COSMO.get_rows(b0, row_range)
  if row_range_col != nothing
    for k in row_range_col
      b_I[k] = b0.nzind[k] + offset
    end
  end

  # also create new set once
  C_new[set_ind] = C
  cone_map[set_ind] = k

  for col = 1:n
    ptr = copy(row_ptr)
    # indices which stores the rows in column for C in A0
    row_range_col = COSMO.get_rows(A0, col, row_range)
    for k in row_range_col
      A_I[k] = A0.rowval[k] +offset
    end
  end


  return row_ptr + C.dim, overlap_ptr, set_ind + 1, sp_ind
end

function clique_rows_map(row_start::Int64, sntree::SuperNodeTree, C::DecomposableCones{<:Real})
  Nc = num_cliques(sntree)
  rows = Array{UnitRange{Int64}}(undef,  Nc)
  ind = zeros(Int64, Nc)
  for iii = Nc:-1:1
    num_rows = COSMO.get_blk_rows(COSMO.get_nBlk(sntree, iii), C)
    rows[iii] = row_start:row_start+num_rows-1
    ind[iii] = sntree.snd_post[iii]
    row_start += num_rows
  end
  return Dict(ind .=> rows)
end

function parent_block_indices(par_clique::Array{Int64, 1}, i::Int64, j::Int64)
  ir = searchsortedfirst(par_clique, i)
  jr = searchsortedfirst(par_clique, j)
  return COSMO.mat_to_svec_ind(ir, jr)
end




function add_entries!(A_I::Array{Int64, 1}, b_I::Array{Int64, 1}, C_new::Array{COSMO.AbstractConvexSet{Float64}, 1}, row_ptr::Int64, A0::SparseMatrixCSC, b0::SparseVector, row_range::UnitRange{Int64}, overlap_ptr::Int64, set_ind::Int64, sp_ind::Int64,
  sp_arr::Array{SparsityPattern, 1},  C::DecomposableCones{<: Real}, k::Int64,  cone_map::Dict{Int64, Int64})

  sp = sp_arr[sp_ind]
  sntree = sp.sntree
  ordering = sp.ordering

  m, n = size(A0)

  # determine the row ranges for each of the subblocks
  clique_to_rows = COSMO.clique_rows_map(row_ptr, sntree, C)

  # loop over cliques in descending topological order

  for iii = num_cliques(sntree):-1:1

    sep = map(v -> ordering[v], get_sep(sntree, iii))
    isa(sep, Array{Any, 1}) && (sep = Int64[])
    clique = map(v -> ordering[v], get_clique(sntree, iii))
    sort!(clique)

    if iii == num_cliques(sntree)
      par_clique = Int64[]
      par_rows = 0:0
    else
      par_ind = COSMO.get_clique_par(sntree, iii)
      par_rows = clique_to_rows[par_ind]
      par_clique = map(v -> ordering[v], get_clique_by_ind(sntree, par_ind))
      sort!(par_clique)
    end

    for col = 1:n
      row_range_col = COSMO.get_rows(A0, col, row_range)
      if col == 1
        row_range_b = COSMO.get_rows(b0, row_range)
      else
        row_range_b = 0:0
      end
      overlap_ptr = add_clique_entries!(A_I, b_I, A0, b0, clique, sep, par_clique, par_rows, col, C, row_ptr, overlap_ptr, row_range, row_range_col, row_range_b)
    end

    num_rows = get_blk_rows(get_nBlk(sntree, iii), C)
    # create and add new cone for subblock
    cone_map[set_ind] = k
    C_new[set_ind] = typeof(C)(num_rows, sp_ind, iii)
    row_ptr += num_rows
    set_ind += 1
  end
  return row_ptr, overlap_ptr, set_ind, sp_ind + 1
end

function add_clique_entries!(A_I::Array{Int64, 1}, b_I::Array{Int64, 1}, A0::SparseMatrixCSC, b0::SparseVector, clique::Array{Int64, 1}, sep::Array{Int64, 1},  par_clique::Array{Int64, 1}, par_rows::UnitRange{Int64}, col::Int64,  C::DecomposableCones{<: Real}, row_ptr::Int64, overlap_ptr::Int64, row_range::UnitRange{Int64}, row_range_col::UnitRange{Int64}, row_range_b)
  counter = 0
  for j in clique, i in clique
    if i <= j
      if i ∈ sep && j ∈ sep
        if col == 1
          A_I[overlap_ptr] = row_ptr + counter # this creates the +1 entry
          A_I[overlap_ptr + 1] = par_rows.start + COSMO.parent_block_indices(par_clique, i, j) - 1 # this creates the -1 entry
          overlap_ptr += 2
        end
      else
        # find row where i, j is stored in A
        k = COSMO.vectorized_ind(i, j, C.sqrt_dim, C)

        row_0 = COSMO.get_row_index(k, A0.rowval, C.sqrt_dim, row_range, row_range_col)
        # row_0 happens when (i, j) references an edge that was added by merging cliques, the corresponding value will be zero
        # and can be disregarded
        #show(i, j, row_0, row_ptr + counter)
        if row_0 != 0
          A_I[row_0] =  row_ptr + counter
        end

        if col == 1
          row_0 = COSMO.get_row_index(k, b0.nzind, C.sqrt_dim, row_range, row_range_b)
          if row_0 != 0
            b_I[row_0] = row_ptr + counter
          end
        end
      end
    counter += 1
    end
  end
  return overlap_ptr
end


"Given the svec index `k` and an offset `row_range_col.start`, return the location of the (i, j)th entry in the sparse structure A."
function get_row_index(k::Int64, rowval::Array{Int64}, sqrt_dim::Int64, row_range::UnitRange{Int64}, row_range_col::UnitRange{Int64})

  k_shift = row_range.start + k - 1

  # determine upper set boundary of where the row could be
  u = min(row_range_col.stop, row_range_col.start + k_shift - 1)

  # find index of first entry >= k, starting in the interval [l, u]
  # if no, entry is >= k, returns u + 1
  r = searchsortedfirst(rowval, k_shift, row_range_col.start, u, Base.Order.Forward)
  # if no r s.t. rowval[r] = k_shift was found that means that the (i, j)th entry represents an edded edge (zero) from clique merging
  if r > u || rowval[r] != k_shift
    return 0
  else
    return r
  end
end

function alternating_sequence(total_length::Int64, n_start::Int64)
  v = ones(Float64, total_length)
  for i=n_start + 1:2:length(v)
    v[i] = -1
  end
  return v
end

function extra_columns(total_length::Int64, n_start::Int64, start_val::Int64)
  v = zeros(Int64, total_length)
  for i = n_start:2:length(v)-1
    v[i] = start_val
    v[i + 1] = start_val
    start_val += 1
  end
  return v
end

function findnz!(I::Vector{Ti}, J::Vector{Ti}, V::Vector{Tv}, S::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    numnz = nnz(S)
    count = 1
    @inbounds for col = 1 : S.n, k = S.colptr[col] : (S.colptr[col+1]-1)
       # I[count] = S.rowval[k]
        J[count] = col
        V[count] = S.nzval[k]
        count += 1
    end
end

function augment_clique_based!(ws)
  A = ws.p.A
  b = ws.p.b
  q = ws.p.q
  P = ws.p.P
  cones = ws.p.C.sets
  sp_arr = ws.ci.sp_arr


  mA, nA, num_overlapping_entries = find_A_dimension(ws.p.model_size[2], ws.p.C.sets, ws.ci.sp_arr)

  # find number of decomposed and total sets and allocate structure for new compositve convex set
  num_total, num_new_psd_cones = COSMO.num_cone_decomposition(ws)

  C_new = Array{COSMO.AbstractConvexSet{Float64}}(undef, num_total - 1)

  nz = nnz(A)
  Aa_I = zeros(Int64, nz + 2 * num_overlapping_entries)
  Aa_J = extra_columns(nz + 2 * num_overlapping_entries, nz + 1, ws.p.model_size[2] + 1)
  Aa_V =  alternating_sequence(nz + 2 * num_overlapping_entries, nz + 1)
  findnz!(Aa_I, Aa_J, Aa_V, A)

  row_ranges = COSMO.get_set_indices(cones) # the row ranges of each cone in the original problem
  row_ptr = 1 # a continuously updated pointer to the start of the first entry of a cone in A_I
  overlap_ptr = nz + 1 # a continuously updated pointer to the next free entry to enter the rows for the +1, -1 overlap entries

  bs = sparse(b)
  b_I = zeros(Int64, length(bs.nzval))
  b_V = bs.nzval
  sp_ind = 1
  set_ind = 1

  for (k, C) in enumerate(cones)
    row_range = row_ranges[k]
    row_ptr, overlap_ptr, set_ind, sp_ind = COSMO.add_entries!(Aa_I, b_I, C_new, row_ptr, A, bs, row_range, overlap_ptr, set_ind, sp_ind, sp_arr, C, k, ws.ci.cone_map)
  end


  Aa = sparse(Aa_I, Aa_J, Aa_V, mA, nA)
  ba = Vector(SparseArrays._sparsevector!(b_I, b_V, mA))
  #ba = Vector(sparsevec(b_I, b_V, mA))

  ws.p.P = blockdiag(P, spzeros(num_overlapping_entries, num_overlapping_entries))
  ws.p.q = vec([q; zeros(num_overlapping_entries)])
  ws.p.A = Aa
  ws.p.b = ba
  ws.p.model_size[1] = size(ws.p.A, 1)
  ws.p.model_size[2] = size(ws.p.A, 2)
  ws.p.C = COSMO.CompositeConvexSet(C_new)

  return nothing
end

function add_sub_blocks!(s::SplitVector, s_decomp::SplitVector, μ::AbstractVector, μ_decomp::AbstractVector, ci::ChordalInfo, C::CompositeConvexSet, C0::CompositeConvexSet, cone_map::Dict{Int64, Int64})
  sp_arr = ci.sp_arr
  row_start = 1 # the row pointer in the decomposed problem
  row_ranges = get_set_indices(C0.sets) # the row ranges of the same cone (or "parent" cone) in the original problem

  # iterate over all the cones of the decomposed problems and add the entries into the correct positions of the original problem
  for (k, C) in enumerate(C.sets)
    row_range = row_ranges[cone_map[k]]
    row_start = add_blocks!(s, μ, row_start, row_range, sp_arr, s_decomp, μ_decomp, C)
  end
  return nothing
end

function add_blocks!(s::SplitVector, μ::AbstractVector, row_start::Int64, row_range::UnitRange{Int64}, sp_arr::Array{SparsityPattern, 1}, s_decomp::SplitVector, μ_decomp::AbstractVector, C::AbstractConvexSet)

  @. s.data[row_range] = s_decomp.data[row_start:row_start + C.dim - 1]
  @. μ[row_range] = μ_decomp[row_start:row_start + C.dim - 1]
  return row_start + C.dim
end

function add_blocks!(s::SplitVector, μ::AbstractVector, row_start::Int64, row_range::UnitRange{Int64}, sp_arr::Array{SparsityPattern, 1}, s_decomp::SplitVector, μ_decomp::AbstractVector, C::DecomposableCones{ <: Real})
  # load the appropriate sparsity_pattern
  sp = sp_arr[C.tree_ind]
  sntree = sp.sntree
  ordering = sp.ordering
  #row_start = sp.row_range.start
  N = length(ordering)
  clique = map(v -> ordering[v], get_clique(sntree, C.clique_ind))
  sort!(clique)
  counter = 0
  for j in clique, i in clique
    if isa(C, PsdCone) || i <= j
      offset = COSMO.vectorized_ind(i, j, N, C) - 1
      s.data[row_range.start + offset] += s_decomp.data[row_start + counter]
      μ[row_range.start + offset] += μ_decomp[row_start + counter]
      counter += 1
    end

  end
  row_start += get_blk_rows(length(clique), C)

  return row_start
end

function reverse_decomposition!(ws::COSMO.Workspace, settings::COSMO.Settings)

  mO = ws.ci.originalM
  nO = ws.ci.originalN

  vars = Variables{Float64}(mO, nO, ws.ci.originalC)
  vars.x .= ws.vars.x[1:nO]

  if settings.colo_transformation
    # reassemble the original variables s and μ
    add_sub_blocks!(vars.s, ws.vars.s, vars.μ, ws.vars.μ, ws.ci, ws.p.C, ws.ci.originalC, ws.ci.cone_map)
  else
    H = ws.ci.H
    vars.s  .= SplitVector{Float64}(H * ws.vars.s[mO + 1:end], ws.ci.originalC)
    # this performs the operation μ = sum H_k^T *  μ_k which is negative of the (uncompleted) dual variable of the original problem
    vars.μ .= H * ws.vars.μ[mO + 1:end]
  end

  ws.p.C = ws.ci.originalC
  # if user requests, perform positive semidefinite completion on entries of μ that were not in the decomposed blocks
  ws.vars = vars
  settings.complete_dual && psd_completion!(ws)

  return nothing
end

# The psd entries of μ that correspond to the zeros in s are not constrained by the problem
# however, in order to make the dual psd cone positive semidefinite we have to do a
# positive semidefinite completion routine to choose the values
function psd_completion!(ws::COSMO.Workspace)

  # loop over psd cones
  row_ranges = get_set_indices(ws.p.C.sets)
  sp_ind = 1
  for (kkk, C) in enumerate(ws.p.C.sets)
    sp_ind = complete!(ws.vars.μ, C, ws.ci.sp_arr, sp_ind, row_ranges[kkk])
  end

  return nothing
end

complete!(μ::AbstractVector, ::AbstractConvexSet, sp_arr::Array{SparsityPattern}, sp_ind::Int64, rows::UnitRange{Int64}) = sp_ind

function complete!(μ::AbstractVector, C::PsdCone{<: Real}, sp_arr::Array{SparsityPattern}, sp_ind::Int64, rows::UnitRange{Int64})
  sp = sp_arr[sp_ind]

  μ_view = view(μ, rows)
  # make this y = -μ
  @. μ_view *= -1

  M = reshape(μ_view, C.sqrt_dim, C.sqrt_dim)

  psd_complete!(M, C.sqrt_dim, sp.sntree, sp.ordering)

  @. μ_view *= -1
  return sp_ind + 1
end

function complete!(μ::AbstractVector, C::PsdConeTriangle{<: Real}, sp_arr::Array{SparsityPattern}, sp_ind::Int64, rows::UnitRange{Int64})
  sp = sp_arr[sp_ind]

  μ_view = view(μ, rows)

  # I want to psd complete y, which is -μ
  populate_upper_triangle!(C.X, -μ_view, 1. / sqrt(2))
  psd_complete!(C.X, C.sqrt_dim, sp.sntree, sp.ordering)
  extract_upper_triangle!(C.X, μ_view, sqrt(2))
  @. μ_view *= -1
  return sp_ind + 1
end


# positive semidefinite completion (from Vandenberghe - Chordal Graphs..., p. 362)
# input: A - positive definite completable matrix
function psd_complete!(A::AbstractMatrix, N::Int64, sntree::SuperNodeTree, p::Array{Int64})

  # if a clique graph based merge strategy was used for this sparsity pattern, recompute a valid clique tree
  #recompute_clique_tree(sntree.strategy) && clique_tree_from_graph!(sntree, sntree.strategy)

  ip = invperm(p)

  As = Symmetric(A, :U)
  # permutate matrix based on ordering p (p must be a vector type), W is in the order that the cliques are based on
  W = copy(As[p, p])
  W = Matrix(W)
  num_cliques = COSMO.num_cliques(sntree)

  # go through supernode tree in descending order (given a post-ordering). This is ensured in the get_snd, get_sep functions
  for j = (num_cliques - 1):-1:1

    # in order to obtain ν, α the vertex numbers of the supernode are mapped to the new position of the permuted matrix
    # index set of snd(i) sorted using the numerical ordering i,i+1,...i+ni
    ν = get_snd(sntree, j)
    #clique = get_clique(sntree, snd_id)
    # index set containing the elements of col(i) \ snd(i) sorted using numerical ordering σ(i)
    α = get_sep(sntree, j)

    # index set containing the row indizes of the lower-triangular zeros in column i (i: representative index) sorted by σ(i)
    i = ν[1]
    η = collect(i+1:1:N)
    # filter out elements in lower triangular part of column i that are non-zero
    filter!(x -> !in(x, α) && !in(x, ν), η)

    Waa = W[α, α]
    Wαν = view(W, α, ν)
    Wηα = view(W, η, α)

    Y = zeros(length(α), length(ν))
    try
     Y[:, :] = Waa \ Wαν
   catch
    Waa_pinv = pinv(Waa)
    Y[:, :] = Waa_pinv * Wαν
  end

  W[η, ν] =  Wηα * Y
  # symmetry condition
  W[ν, η] = view(W, η, ν)'
end

# invert the permutation
A[:, :] =  W[ip, ip]
end


