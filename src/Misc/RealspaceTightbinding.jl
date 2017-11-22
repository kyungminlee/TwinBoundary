
type RealSpaceTightbindingModel{D, S <: Number}
  orbitals :: Vector{String}
  hoppings :: Vector{HoppingElement{D, S}}

  function TightbindingModel(orbitals :: Vector{String})
    new(orbitals, HoppingElement{D, S}[])
  end

  function TightbindingModel()
    new([:a], HoppingElement{D, S}[])
  end
end
