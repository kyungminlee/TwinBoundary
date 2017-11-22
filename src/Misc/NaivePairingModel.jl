module NaivePairingModel

function make_pairing_momentumspace{T<:Number}(pairing_orderparameter ::T)
  (kx ::Float64, ky ::Float64) -> begin
    ϕ = pairing_orderparameter
    Δ = 4.0 * real(ϕ) * cos(kx/2) * cos(ky/2)
    return Complex128[0 Δ; conj(Δ) 0]
  end
end


function make_pairing_mixedspace{T <: Number}(
      pairing_orderparameter ::T,
      nx ::Int;
      periodic ::Bool=true)
  Φ = copy(pairing_orderparameter)

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
        pairing_gap_matrix[ro, rs, co, cs] += v
        pairing_gap_matrix[co, cs, ro, rs] += conj(v)
      end
    end

    collect = periodic? collect_periodic : collect_open

    emy = exp(-0.5im * ky)
    epy = exp( 0.5im * ky)
    for ix=1:nx
      # A-type plaquette
      let ϕ = pairing_orderparameter,
          vmy = 0.5 * ϕ * emy,
          vpy = 0.5 * ϕ * epy
        collect(1, ix  , 2, ix  , vpy)
        collect(2, ix  , 1, ix+1, vmy)
        collect(1, ix+1, 2, ix  , vmy)
        collect(2, ix  , 1, ix  , vpy)
      end

      # B-type plaquette
      let ϕ = pairing_orderparameter,
          vmy = 0.5 * ϕ * emy,
          vpy = 0.5 * ϕ * epy
        collect(2, ix  , 1, ix+1, vpy)
        collect(1, ix+1, 2, ix+1, vmy)
        collect(2, ix+1, 1, ix+1, vmy)
        collect(1, ix+1, 2, ix  , vpy)
      end
    end
    return reshape(pairing_gap_matrix, (2*nx, 2*nx))
  end
end


function make_pairing_realspace_dense(pairing_orderparameter ::Number,
                                      nx ::Int, ny ::Int;
                                      periodic ::Bool=true)
  pairing_gap_matrix = zeros(Complex128, (2, nx, ny, 2, nx, ny))

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
    v = 0.5 * pairing_orderparameter
    collect(1, ix  , iy  , 2, ix  , iy  , v)
    collect(2, ix  , iy  , 1, ix+1, iy  , v)
    collect(1, ix+1, iy  , 2, ix  , iy-1, v)
    collect(2, ix  , iy-1, 1, ix  , iy  , v)

    collect(1, ix  , iy  , 2, ix-1, iy  , v)
    collect(2, ix-1, iy  , 1, ix  , iy+1, v)
    collect(1, ix  , iy+1, 2, ix  , iy  , v)
    collect(2, ix  , iy  , 1, ix  , iy  , v)
  end
  return reshape(pairing_gap_matrix, (2*nx*ny, 2*nx*ny))
end


end