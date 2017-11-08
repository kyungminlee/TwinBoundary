using Kore

#=
A: site A
B: site B
1: Plaquette 1
2: Plaquette 2

1 B 1 B 1 B 1 B
A 2 A 2 A 2 A 2
1 B 1 B 1 B 1 B
A 2 A 2 A 2 A 2
1 B 1 B 1 B 1 B
A 2 A 2 A 2 A 2
1 B 1 B 1 B 1 B
=#

function make_pairing_momentumspace{T<:Number}(
      pairing_orderparameter ::T,
      formfactor)
  const R = [Float64[0.0, 0.0], Float64[0.5, 0.5]]
  (kx ::Float64, ky ::Float64) -> begin
    ret = zeros(Complex128, (2, 2))
    for i_plq in [1,2]
      for (ro, rx, ry, co, cx, cy, vff) in formfactor[i_plq]
        dis = [cx, cy] - [rx, ry] + R[co] - R[ro]
        v = 0.5 * vff * pairing_orderparameter * exp(1im * dis ⋅ [kx, ky])
        ret[ro, co] += v
        ret[co, ro] += conj(v)
      end
    end
    return ret
  end
end

# 1 : nearest neighbor

const FORMFACTORS = Dict{String, Any}(
  "s1" => Dict(
    1 => [
      ( 1, 0, 0, 2,-1, 0, 1),
      ( 2,-1, 0, 1, 0, 1, 1),
      ( 1, 0, 1, 2, 0, 0, 1),
      ( 2, 0, 0, 1, 0, 0, 1)
    ],

    2 => [
      (1, 0, 0, 2, 0,  0, 1),
      (2, 0, 0, 1, 1,  0, 1),
      (1, 1, 0, 2, 0, -1, 1),
      (2, 0,-1, 1, 0,  0, 1),
    ],
  ),

  "d1" => Dict(
    1 => [
      ( 1, 0, 0, 2,-1, 0, -1),
      ( 2,-1, 0, 1, 0, 1,  1),
      ( 1, 0, 1, 2, 0, 0, -1),
      ( 2, 0, 0, 1, 0, 0,  1)
    ],
    2 => [
      (1, 0, 0, 2, 0,  0,  1),
      (2, 0, 0, 1, 1,  0, -1),
      (1, 1, 0, 2, 0, -1,  1),
      (2, 0,-1, 1, 0,  0, -1),
    ],
  ),
)


function make_pairing_realspace_dense{T <: Number}(
      pairing_orderparameter :: Array{T, 3},
      formfactor;
      periodic ::Bool=true)

  (n_basis, nx, ny) = size(pairing_orderparameter)
  @assert n_basis == 2

  pairing_gap_matrix = zeros(T, (2, nx, ny, 2, nx, ny))

  collect_periodic = (ro, rx, ry, co, cx, cy, v) -> begin
    rx, ry = mod(rx, nx, 1), mod(ry, ny, 1)
    cx, cy = mod(cx, nx, 1), mod(cy, ny, 1)
    pairing_gap_matrix[ro, rx, ry, co, cx, cy] += v
    pairing_gap_matrix[co, cx, cy, ro, rx, ry] += conj(v)
  end

  collect_open = (ro, rx, ry, co, cx, cy, v) -> begin
    if (1 <= rx <= nx && 1 <= cx <= nx &&
        1 <= ry <= ny && 1 <= cy <= ny)
      pairing_gap_matrix[ro, rx, ry, co, cx, cy] += v
      pairing_gap_matrix[co, cx, cy, ro, rx, ry] += conj(v)
    end
  end

  collect = periodic ? collect_periodic : collect_open

  for ix=1:nx, iy=1:ny
    for i_plq in [1,2]
      v = 0.5 * pairing_orderparameter[i_plq, ix, iy]
      for (ro, rx, ry, co, cx, cy, vff) in formfactor[i_plq]
        collect(ro, ix+rx, iy+ry, co, ix+cx, iy+cy, vff*v)
      end
    end
  end
  return reshape(pairing_gap_matrix, (2*nx*ny, 2*nx*ny))
end


function make_pairing_mixedspace{T <: Number}(
      pairing_orderparameter ::Array{T, 2},
      formfactor;
      periodic ::Bool=true)

  (n_basis, nx) = size(pairing_orderparameter)
  @assert(n_basis == 2, "number of orbitals should be 2")

  Φ = copy(pairing_orderparameter)

  const orbital_displacement = [0.0, 0.5]

  (ky ::Float64) -> begin
    # two for basis in doubled unit cell
    pairing_gap_matrix = zeros(Complex128, (2, nx, 2, nx))

    collect_periodic = (ro, ru, co, cu, v) -> begin
      ru, cu = mod(ru, nx, 1), mod(cu, nx, 1)
      pairing_gap_matrix[ro, ru, co, cu] += v
      pairing_gap_matrix[co, cu, ro, ru] += conj(v)
    end

    collect_open = (ro, ru, co, cu, v) -> begin
      if 1 <= ru <= nx && 1 <= cu <= nx
        pairing_gap_matrix[ro, ru, co, cu] += v
        pairing_gap_matrix[co, cu, ro, ru] += conj(v)
      end
    end

    collect = periodic? collect_periodic : collect_open

    for i_plq in [1,2]
      for (ro, rx, ry, co, cx, cy, vff) in formfactor[i_plq]
        dy = (cy - ry) + orbital_displacement[co] - orbital_displacement[ro]

        vff_ky = exp(1im * dy * ky) * vff
        for ix=1:nx
          v = 0.5 * vff_ky * pairing_orderparameter[i_plq, ix]
          collect(ro, ix+rx, co, ix+cx, v)
        end
      end
    end
    return reshape(pairing_gap_matrix, (2*nx, 2*nx))
  end
end
