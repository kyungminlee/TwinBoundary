
function foo()
  plaquette_corners = [(0,0), (1,0), (1,1), (0,1)]

  bonds = Dict()

  # Open boundary condition
  collect_open = (ix, iy, jx, jy, v) -> begin
      half_nx = div(nx, 2)
      half_ny = div(ny, 2)
      if (1 <= ix <= nx && 1 <= iy <= ny &&
          1 <= jx <= nx && 1 <= jy <= ny)
          dx = mod(jx - ix + half_nx, nx) - half_nx
          dy = mod(jy - iy + half_ny, ny) - half_ny
          bonds[(ix, iy, jx, jy, dx, dy)] = get(bonds, (ix, iy, jx, jy, dx, dy), 0.0 + 0.0im) + v
      end
  end

  # Periodic
  collect_periodic = (ix, iy, jx, jy, v) -> begin
      half_nx = div(nx, 2)
      half_ny = div(ny, 2)
      ix, iy = mod(ix, nx, 1), mod(iy, ny, 1)
      jx, jy = mod(jx, nx, 1), mod(jy, ny, 1)
      dx = mod(jx - ix + half_nx, nx) - half_nx
      dy = mod(jy - iy + half_ny, ny) - half_ny
      bonds[(ix, iy, jx, jy, dx, dy)] = get(bonds, (ix, iy, jx, jy, dx, dy), 0.0 + 0.0im) + v
  end

  #collect = collect_periodic
  collect = collect_open

  for (dx, dy, vff) in [( 1, 0, Δs),(-1, 0, Δs),( 0, 1, Δs),( 0,-1, Δs)]
    map_bonds = [(0,0), (1,0), (1,1), (0,1)]
    for px=1:nx, py=1:ny
      dx1, dy1 = map_bonds[1]
      dx2, dy2 = map_bonds[2]
      dx3, dy3 = map_bonds[3]
      dx4, dy4 = map_bonds[4]

      ix, iy = (px + dx1), (py + dy1)
      jx, jy = (px + dx2), (py + dy2)
      kx, ky = (px + dx3), (py + dy3)
      lx, ly = (px + dx4), (py + dy4)

      v = 0.25 * vff * map_orderparameter[px, py]
      # TODO: HERE!!!

      collect(ix, iy, jx, jy, v)
      collect(jx, jy, kx, ky, v)
      collect(kx, ky, lx, ly, v)
      collect(lx, ly, ix, iy, v)

      collect(ix, iy, lx, ly, conj(v))
      collect(lx, ly, kx, ky, conj(v))
      collect(kx, ky, jx, jy, conj(v))
      collect(jx, jy, ix, iy, conj(v))
    end
  end

  for (dx, dy, vff) in [( 1, 0, Δd),(-1, 0, Δd),( 0, 1,-Δd),( 0,-1,-Δd)]
    for ix=1:nx, iy=1:ny
      jx = ix + dx
      jy = iy + dy
      collect(ix, iy, jx, jy, vff)
    end
  end

end
