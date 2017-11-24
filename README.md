

```julia
type TightbindingModel{D, S <: Number}
    orbitals :: Vector{String}
    hoppings :: Vector{HoppingElement{D, S}}

    ...
end

function add_hopping!{D, S<:Number, S2 <: Number}(
  tightbinding :: TightbindingModel{D, S},
  displacement :: NTuple{D, Int},
  row_orbital :: String,
  col_orbital :: String,
  value :: S2)

function conjugate{D, S <: Number}(tightbinding :: TightbindingModel{D, S})

function isHermitian{D, S}(tb_model :: TightbindingModel{D, S}, epsilon ::Real = 1E-10)

function doubleUnitcell{S}(tb_model ::TightbindingModel{2, S})

function makeMixedspace{S}(tb_model ::TightbindingModel{2, S},
                           ny       ::Real,
                           periodic ::Bool = false)
function makeMixedspace2{D, S}(tb_model ::TightbindingModel{D, S},
                               size :: NTuple{D, Int})
function makeRealspaceDense{D, S}(tb_model::TightbindingModel{D, S},
                                  size ::NTuple{D, Int},
                                  periodic ::Bool = true)
function makeRealspaceSparse{D, S}(tb_model::TightbindingModel{D, S},
                                   size ::NTuple{D, Int},
                                   periodic ::Bool = true)

```