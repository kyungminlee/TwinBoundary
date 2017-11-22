using Base.Test

include("Tightbinding.jl")
using Tightbinding

Base.ishermitian{T <: Number}(m ::Matrix{T}; atol=0.0) = isapprox(m, m'; atol=atol)

#let
  t = 1.0 + 0.5im
  tb = TightbindingModel{2, Complex128}([("s", [0.0, 0.0])])

  @show tb

  add_hopping!(tb, ( 1, 0), "s", "s", t)
  add_hopping!(tb, ( 0, 1), "s", "s", t)

  @test !ishermitian(tb)

  add_hopping!(tb, (-1, 0), "s", "s", conj(t))
  add_hopping!(tb, ( 0,-1), "s", "s", conj(t))

  @test ishermitian(tb)

  nx, ny = 3, 5

  kxs = linspace(0, 2 * pi, nx+1)[1:(end-1)]
  kys = linspace(0, 2 * pi, ny+1)[1:(end-1)]  

  h_real_dense = make_realspace_dense(tb, (nx, ny))
  h_real_sparse = make_realspace_sparse(tb, (nx, ny))
  h_momentum = make_momentumspace(tb)
  h_mixed1 = make_mixedspace(tb, nx; periodic=true)
  h_mixed2 = make_mixedspace2(tb, (nx, 0))

  @test h_real_dense ≈ full(h_real_sparse)
  
  e_real_dense = eigvals( Hermitian(0.5 * (h_real_dense + h_real_dense')) )
  e_momentum = Float64[]
  e_mixed = Float64[]
  for iy=1:ny
    ky = kys[iy]
    for ix=1:nx
      kx = kxs[ix]
      hk = h_momentum(kx, ky)
      @test ishermitian(hk; atol=1E-8)
      append!(e_momentum, eigvals(Hermitian(0.5 * (hk + hk'))))
    end
    hm1 = h_mixed1(ky)
    hm2 = h_mixed2([ky])
    @test hm1 ≈ hm2
    @test ishermitian(hm1)
    append!(e_mixed, eigvals(hm1))
  end
  sort!(e_real_dense)
  sort!(e_momentum)
  sort!(e_mixed)

  @test e_momentum ≈ e_real_dense
  @test e_momentum ≈ e_mixed

  tbc = conjugate(tb)
  @test ishermitian(tbc)
  
  h_real_dense_conj = make_realspace_dense(tbc, (nx, ny))
  @test ishermitian(h_real_dense_conj)
  e_real_dense_conj = eigvals(Hermitian(0.5*(h_real_dense_conj + h_real_dense_conj')))
  sort!(e_real_dense_conj)
  @test e_momentum ≈ e_real_dense_conj

#end