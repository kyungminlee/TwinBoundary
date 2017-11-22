using Base.Test
using PyPlot
using Kore

include("Tightbinding.jl")
include("PairingModel.jl")

nx, ny = 4, 4
kxs = linspace(0, 2*pi, nx+1)[1:end-1]
kys = linspace(0, 2*pi, ny+1)[1:end-1]

#pairing_orderparameter = ones(Complex128, (2, nx))
Δ = 1.0 + 0.3im

pairing_realspace_dense = make_pairing_realspace_dense(Δ, nx, ny)
pairing_momentumspace = make_pairing_momentumspace(Δ)

@show maximum(abs(pairing_realspace_dense - pairing_realspace_dense'))

@test isapprox(pairing_realspace_dense, pairing_realspace_dense'; atol=1E-8)

energy_realspace = eigvals(Hermitian(pairing_realspace_dense))
energy_momentumspace = Float64[]
for kx in kxs, ky in kys
  Δk = pairing_momentumspace(kx, ky)
  @test ishermitian(Δk)
  append!(energy_momentumspace, eigvals(Δk))
end

sort!(energy_realspace)
sort!(energy_momentumspace)

chopzero!(energy_realspace; tol=1E-8)
chopzero!(energy_momentumspace; tol=1E-8)

#@show energy_realspace
#@show energy_momentumspace
@show isapprox(energy_realspace, energy_momentumspace)


one_dimensional_orderparameter = Δ * ones(Complex128, (2, nx))
pairing_mixedspace = make_pairing_mixedspace(one_dimensional_orderparameter)

@test all([ishermitian(pairing_mixedspace(ky)) for ky in kys])

energy_mixedspace = vcat([eigvals(pairing_mixedspace(ky)) for ky in kys]...)
sort!(energy_mixedspace)
chopzero!(energy_mixedspace)
#@show energy_mixedspace
@show isapprox(energy_realspace, energy_mixedspace)


#=
Comparing momentum space version with real space version.
Note that the real space fourier transform is carried over unit cells, not sites. (i.e. zero gauge connection between two basis vectors within a unit cell.), and thus the blocks can differ by a U(1) gauge transformation.
=#
let pp = reshape(pairing_realspace_dense, (2, nx, ny, 2, nx, ny))
  pp2 = fft( ifft(pp, (2, 3)), (5, 6))

  datas = []
  for kx in kxs, ky in kys
    push!(datas, pairing_momentumspace(kx, ky))
  end
  using JLD
  save("compare.jld",
        "real", reshape(pp2, (2*nx*ny, 2*nx*ny)),
        "momentum", datas)
end