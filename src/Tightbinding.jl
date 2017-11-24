module Tightbinding

using Kore

export TightbindingModel
export add_hopping!
export conjugate
export double_unitcell
export ishermitian
export make_realspace_dense
export make_realspace_sparse
export make_momentumspace
export make_mixedspace
export make_mixedspace2
export nambufy


immutable HoppingElement{D, S <: Number}
  displacement :: NTuple{D, Int}
  row_orbital  :: String
  col_orbital  :: String
  value        :: S
end

function Base.show(io ::IO, he ::HoppingElement)
  print(io, "HE(")
  print(io, "d: $(he.displacement), ")
  print(io, "r: $(he.row_orbital), ")
  print(io, "c: $(he.col_orbital), ")
  print(io, "v: $(he.value))")
end

function convert{T <: Tuple, D, S}(::Type{T}, he::HoppingElement{D, S})
  (he.displacement, he.row_orbital, he.col_orbital, he.value)
end


function convert{D, R}(::HoppingElement{D, Complex{R}}, he::HoppingElement{D, R})
  return HoppingElement{D, Complex{R}}(he.displacement, he.row_orbital, he.col_orbital, convert(Complex{R}, he.value))
end

immutable OrbitalInformation
  index :: Int
  position :: Vector{Float64}
end

function Base.show(io::IO, obj ::OrbitalInformation)
  print(io, "OI(index:$(obj.index), position:$(obj.position))")
end


mutable struct TightbindingModel{D, S <: Number}
  orbitals :: Vector{String}
  orbital_info :: Dict{String, OrbitalInformation}
  hoppings :: Vector{HoppingElement{D, S}}

  function TightbindingModel{D, S}(
      orbitals :: Vector{String},
      orbital_info :: Dict{String, OrbitalInformation},
      hoppings :: Vector{HoppingElement{D, S}}) where {D, S<:Number}
    new{D, S}(orbitals, orbital_info, hoppings)
  end
end

function TightbindingModel{D, S}(orbitals ::Vector{Tuple{String, Vector{Float64} } } ) where {D,S}
  orbs = String[]
  orbital_info = Dict{String, OrbitalInformation}()
  for (idx, (orb, pos)) in enumerate(orbitals)
    push!(orbs, orb)
    orbital_info[orb] = OrbitalInformation(idx, pos)
  end
  hoppings = HoppingElement{D, S}[]
  #@show orbs
  #@show orbital_info
  TightbindingModel{D, S}(orbs, orbital_info, hoppings)
end

function TightbindingModel{D, S}() where {D,S}
  TightbindingModel{D,S}([("s", zeros(Float64, D))])
end

function Base.show{D, S <: Number}(io ::IO, tb ::TightbindingModel{D, S})
  println(io, "orbitals : ", tb.orbitals)
  println(io, "orbital_info : ", tb.orbital_info)
  println(io, "hoppings :")
  println(io, tb.hoppings)
end

function add_hopping!{D, S<:Number, S2 <: Number}(
  tightbinding :: TightbindingModel{D, S},
  displacement :: NTuple{D, Int},
  row_orbital :: String,
  col_orbital :: String,
  value :: S2)

  if !(row_orbital in tightbinding.orbitals)
    throw(DomainError())#"orbital $(row_orbital) not in tightbinding orbitals ($(tightbinding.orbitals))"))
  end
  if !(col_orbital in tightbinding.orbitals)
    throw(DomainError())#"orbital $(col_orbital) not in tightbinding orbitals ($(tightbinding.orbitals))"))
  end
  if length(displacement) != D
    throw(DomainError())#"spatial displacement should be $D dimensional"))
  end

  push!(tightbinding.hoppings,
        HoppingElement(displacement,
                       row_orbital, col_orbital,
                       Base.convert(S, value)) )
end


function conjugate{D, S <: Number}(tightbinding :: TightbindingModel{D, S})
  if S <: Real
    return tightbinding
  else
    hoppings = HoppingElement{D, S}[
      HoppingElement(hop.displacement,
                     hop.row_orbital,
                     hop.col_orbital,
                     conj(hop.value))
        for hop in tightbinding.hoppings
    ]
    new_tb = TightbindingModel{D, S}(tightbinding.orbitals, tightbinding.orbital_info, hoppings)
    return new_tb
  end
end


"""
old orbitals : [o1, o2, o3]
new orbitals : [o1-a, o2-a, o3-a, o1-b, o2-b, o3-b]
"""
function double_unitcell{S}(tb_model ::TightbindingModel{2, S})
  new_orbitals = Tuple{String, Vector{Float64}}[]
  for orb in tb_model.orbitals
    pos = tb_model.orbital_info[orb].position
    new_pos = [0.5 -0.5; 0.5 0.5] * pos
    push!(new_orbitals, ("$orb-a", new_pos))
  end

  for orb in tb_model.orbitals
    pos = tb_model.orbital_info[orb].position
    new_pos = [0.5 -0.5; 0.5 0.5] * pos
    push!(new_orbitals, ("$orb-b", new_pos + [0.5, 0.5]))
  end

  new_tb_model = TightbindingModel{2, S}(new_orbitals)
  for hop in tb_model.hoppings
    (dx, dy) = hop.displacement
    ro = hop.row_orbital
    co = hop.col_orbital
    val = hop.value
    Dx, Dy = (dx-dy), (dx+dy)

    if mod(Dx, 2) == 0
      add_hopping!(new_tb_model, (Dx ÷ 2, Dy ÷ 2), "$ro-a", "$co-a", val)
      add_hopping!(new_tb_model, (Dx ÷ 2, Dy ÷ 2), "$ro-b", "$co-b", val)
    else
      add_hopping!(new_tb_model, ((Dx - 1) ÷ 2, (Dy - 1) ÷ 2), "$ro-a", "$co-b", val)
      add_hopping!(new_tb_model, ((Dx + 1) ÷ 2, (Dy + 1) ÷ 2), "$ro-b", "$co-a", val)
    end
  end
  return new_tb_model
end

function nambufy{S}(tb_model ::TightbindingModel{2, S})
  new_orbitals = Tuple{String, Vector{Float64}}[]
  for orb in tb_model.orbitals
    pos = tb_model.orbital_info[orb].position
    push!(new_orbitals, ("$orb-p", pos))
    push!(new_orbitals, ("$orb-h", pos))
  end

  new_tb_model = TightbindingModel{2, S}(new_orbitals)
  for hop in tb_model.hoppings
    (dx, dy) = hop.displacement
    ro = hop.row_orbital
    co = hop.col_orbital
    val = hop.value
    add_hopping!(new_tb_model, (dx, dy), "$ro-p", "$co-p", val)
    add_hopping!(new_tb_model, (-dx, -dy), "$co-h", "$ro-h", -val)
  end
  return new_tb_model
end

import Base: ishermitian

function Base.ishermitian{D, S}(tb_model ::TightbindingModel{D, S}; epsilon ::Real = 1E-10)
  orbital_info = tb_model.orbital_info
  hoppings = tb_model.hoppings

  diff = Dict{NTuple{D+2, Int}, S}()
  for hop in hoppings
    dis = hop.displacement

    # diagonal element
    if all([d == 0 for d in dis]) && iro == ico
      if !isreal(v)
        return false
      end
    end

    mdis = -[dis...]
    ro, co = hop.row_orbital, hop.col_orbital
    iro, ico = orbital_info[ro].index, orbital_info[co].index

    k1 = (dis..., iro, ico)
    k2 = (mdis..., ico, iro)

    diff[k1] = haskey(diff, k1) ? (diff[k1] + hop.value) : +hop.value
    diff[k2] = haskey(diff, k2) ? (diff[k2] - conj(hop.value)) : -conj(hop.value)
  end

  for (k, v) in diff
    if isapprox(v, zero(D); atol=epsilon)
      delete!(diff, k)
    end
  end

  return isempty(diff)
end


function make_mixedspace{S}(tb_model ::TightbindingModel{2, S},
                            ny       ::Int;
                            periodic ::Bool=true)
  n_orbital = length(tb_model.orbitals)
  hoppings = tb_model.hoppings
  orbital_info = copy(tb_model.orbital_info)

  (kx ::Real) -> begin
    m = zeros(Complex128, n_orbital, ny, n_orbital, ny)
    for hop in hoppings
      dis = hop.displacement
      ro, co, val = hop.row_orbital, hop.col_orbital, hop.value
      iro, ico = orbital_info[ro].index, orbital_info[co].index
      rpos = orbital_info[ro].position
      cpos = orbital_info[co].position

      # TODO: CHECK SIGN OF CO AND RO. WHICH IS RIGHT?
      kphase = kx ⋅ (dis[1] + cpos[1] - rpos[1])
      v = hop.value * exp(1im * kphase)

      for iy in 1:ny
        jy = iy + dis[2]
        if periodic
          jy = mod(jy, ny, 1)
          m[iro, iy, ico, jy] += v
        else
          if 1 <= jy <= ny
            m[iro, iy, ico, jy] += v
          end # if
        end # if periodic
      end # if for
    end # for hop
    return reshape(m, (n_orbital * ny, n_orbital * ny))
  end # function of k
end


# zero if momentum space
function make_mixedspace2{D, S}(tb_model ::TightbindingModel{D, S},
                               size :: NTuple{D, Int})
  n_orbital = length(tb_model.orbitals)
  hoppings = tb_model.hoppings
  orbital_info = copy(tb_model.orbital_info)

  realspace_indices = zeros(Bool, D)
  momentumspace_indices = zeros(Bool, D)

  for (i, s) in enumerate(size)
    if s != 0
      realspace_indices[i] = true
    else
      momentumspace_indices[i] = true
    end
  end

  Dr = sum(realspace_indices)
  Dk = sum(momentumspace_indices)

  (k ::Vector{Float64}) -> begin
    m = zeros(Complex128, (n_orbital, size[realspace_indices]...,
                           n_orbital, size[realspace_indices]...))
    for hop in hoppings
      dis = [hop.displacement...]
      ro, co, val = hop.row_orbital, hop.col_orbital, hop.value
      iro, ico = orbital_info[ro].index, orbital_info[co].index
      rpos = orbital_info[ro].position
      cpos = orbital_info[co].position
      #TODO CHECK SIGN
      kphase = [k...] ⋅ (
        dis[momentumspace_indices]
        + cpos[momentumspace_indices] - rpos[momentumspace_indices])
      v = hop.value * exp(1im * kphase)
      for i in Base.product([1:n for n in size[realspace_indices]]...)
        j = (([i...] + dis[realspace_indices])...)
        j = mod(j, size[realspace_indices], 1)
        m[iro, i..., ico, j...] += v
      end
    end
    n = n_orbital * prod(size[realspace_indices])
    return reshape(m, (n, n))
  end
end

"""
    make_momentumspace{D, S}(tb_model::TightbindingModel{D, S})

Returns a function of \$\\mathbf{k}\$ which returns the Hamiltonian
at momentum \$\\mathbf{k}\$.
"""
function make_momentumspace{D, S}(tb_model::TightbindingModel{D, S})
  n_orbital = length(tb_model.orbitals)
  hoppings = tb_model.hoppings
  orbital_info = copy(tb_model.orbital_info)

  (k ::Vararg{Real, D}) -> begin
    m = zeros(Complex128, n_orbital, n_orbital)
    for hop in hoppings
      ro = hop.row_orbital
      co = hop.col_orbital
      iro, ico = orbital_info[ro].index, orbital_info[co].index
      rpos = orbital_info[ro].position
      cpos = orbital_info[co].position

      kphase = [k...] ⋅ ([hop.displacement...] + cpos - rpos)
      v = hop.value * exp(1im * kphase)
      m[iro, ico] += v
      #if iro != ico
      #  m[ico, iro] += conj(v)
      #end
    end
    return m
  end
end


function make_realspace_dense{D, S}(tb_model::TightbindingModel{D, S},
                                    size ::NTuple{D, Int};
                                    periodic ::Bool = true)
  n_orbital = length(tb_model.orbitals)
  hoppings = tb_model.hoppings
  orbital_info = tb_model.orbital_info

  hamiltonian = zeros(S, n_orbital, size..., n_orbital, size...)
  for hop in hoppings
    dis = hop.displacement
    ro, co, val = hop.row_orbital, hop.col_orbital, hop.value
    iro, ico = orbital_info[ro].index, orbital_info[co].index
    for i in Base.product([1:n for n in size]...)
      j = (([i...] + [dis...])...)
      if periodic
        j = mod(j, size, 1)
        hamiltonian[iro, i..., ico, j...] += val
      else
        if within(j, size)
          hamiltonian[iro, i..., ico, j...] += val
        end
      end
    end
  end
  len = n_orbital * prod(size)
  reshape(hamiltonian, (len, len))
end


function make_realspace_sparse{D, S}(tb_model::TightbindingModel{D, S},
                                     size ::NTuple{D, Int};
                                     periodic ::Bool = true)
  n_orbital = length(tb_model.orbitals)
  hoppings = tb_model.hoppings
  orbital_info = tb_model.orbital_info

  rows, cols, vals = Int[], Int[], S[]
  collect = (iro, irs, ico, ics, val) -> begin
    r = sub2ind((n_orbital, size...), iro, irs...)
    c = sub2ind((n_orbital, size...), ico, ics...)
    push!(rows, r)
    push!(cols, c)
    push!(vals, val)
  end

  for hop in hoppings
    dis = hop.displacement
    ro, co, val = hop.row_orbital, hop.col_orbital, hop.value
    iro, ico = orbital_info[ro].index, orbital_info[co].index
    for i in Base.product([1:n for n in size]...)
      j = (([i...] + [dis...])...)
      if periodic
        j = mod(j, size, 1)
        collect(iro, i, ico, j, val)
      else
        if within(j, size)
          collect(iro, i, ico, j, val)
        end
      end
    end
  end
  len = n_orbital * prod(size)
  return sparse(rows, cols, vals, len, len)
end


end
