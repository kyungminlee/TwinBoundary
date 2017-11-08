module DaghoferModel

using Formatting
using Tightbinding

export DaghoferParameter
#=
export makeMomentumspace
export makeRealspaceDifference
export makeRealspaceDense
export makeRealspaceSparse
=#

type DaghoferParameter
  t1  ::Float64
  t2  ::Float64
  t3  ::Float64
  t4  ::Float64
  t5  ::Float64
  t6  ::Float64
  t7  ::Float64
  t8  ::Float64
  Δxy ::Float64
  μ   ::Float64
end

function hoppingList(param ::DaghoferParameter)
    t1 = param.t1
    t2 = param.t2
    t3 = param.t3
    t4 = param.t4
    t5 = param.t5
    t6 = param.t6
    t7 = param.t7
    t8 = param.t8
    Δxy = param.Δxy
    μ   = param.μ
    hl = [
        (( 1, 0), "xz", "xz", t1),    ((-1, 0), "xz", "xz", t1),
        (( 0, 1), "xz", "xz", t2),    (( 0,-1), "xz", "xz", t2),
        (( 1, 1), "xz", "xz", t3),    (( 1,-1), "xz", "xz", t3),
        ((-1, 1), "xz", "xz", t3),    ((-1,-1), "xz", "xz", t3),
        (( 0, 0), "xz", "xz", -μ),

        (( 1, 0), "yz", "yz", t2),    ((-1, 0), "yz", "yz", t2),
        (( 0, 1), "yz", "yz", t1),    (( 0,-1), "yz", "yz", t1),
        (( 1, 1), "yz", "yz", t3),    (( 1,-1), "yz", "yz", t3),
        ((-1, 1), "yz", "yz", t3),    ((-1,-1), "yz", "yz", t3),
        (( 0, 0), "yz", "yz", -μ),

        (( 1, 0), "xy", "xy", t5),    ((-1, 0), "xy", "xy", t5),
        (( 0, 1), "xy", "xy", t5),    (( 0,-1), "xy", "xy", t5),
        (( 1, 1), "xy", "xy", t6),    (( 1,-1), "xy", "xy", t6),
        ((-1, 1), "xy", "xy", t6),    ((-1,-1), "xy", "xy", t6),
        (( 0, 0), "xy", "xy", Δxy-μ),

        (( 1, 1), "xz", "yz",-t4),    (( 1,-1), "xz", "yz", t4),
        ((-1, 1), "xz", "yz", t4),    ((-1,-1), "xz", "yz",-t4),

        (( 1, 0), "xz", "xy",-t7),    ((-1, 0), "xz", "xy", t7),

        (( 0, 1), "yz", "xy",-t7),    (( 0,-1), "yz", "xy", t7),

        (( 1, 1), "xz", "xy",-t8),    ((-1, 1), "xz", "xy", t8),
        (( 1,-1), "xz", "xy",-t8),    ((-1,-1), "xz", "xy", t8),

        (( 1, 1), "yz", "xy",-t8),    ((-1, 1), "yz", "xy",-t8),
        (( 1,-1), "yz", "xy", t8),    ((-1,-1), "yz", "xy", t8),
    ]
    hl2 = copy(hl)
    for ((dx, dy), ro, co, v) in hl
        if ro != co
            push!(hl2, ((-dx, -dy), co, ro, conj(v)))
        end
    end

    hl2
end

function daghoferModel(S :: Type, param :: DaghoferParameter)
    @assert( S <: Number )
    tb_model = TightbindingModel{2, S}([
        ("xz", [0.0, 0.0]),
        ("yz", [0.0, 0.0]) ,
        ("xy", [0.0, 0.0])])

    hops = hoppingList(param)
    for (dis, ro, co, v) in hops
        add_hopping!(tb_model, dis, ro, co, v)
    end
    return tb_model
end

function spinfulDaghoferModel(S ::Type,
                              param ::DaghoferParameter,
                              soc::Float64)
    @assert( S <: Number )
    tb_model = TightbindingModel{2, S}([
        ("xz_up", [0.0, 0.0]),
        ("yz_up", [0.0, 0.0]),
        ("xy_up", [0.0, 0.0]),
        ("xz_dn", [0.0, 0.0]),
        ("yz_dn", [0.0, 0.0]),
        ("xy_dn", [0.0, 0.0]),
    ])

    # spin-independent hoppings
    hops = hoppingList(param)
    for sp in ["up", "dn"], (dis, ro, co, v) in hops
        ros = format("{}_{}", ro, sp)
        cos = format("{}_{}", co, sp)
        add_hopping!(tb_model, dis, ros, cos, v)
    end
    # spin-orbit coupling
    add_hopping!(tb_model, (0,0), "xz_up", "yz_up", +im * soc)
    add_hopping!(tb_model, (0,0), "yz_up", "xz_up", -im * soc)
    add_hopping!(tb_model, (0,0), "xz_dn", "yz_dn", -im * soc)
    add_hopping!(tb_model, (0,0), "yz_dn", "xz_dn", +im * soc)

    return tb_model
end

#=
function makeMomentumspace(param::DaghoferParameter)
    t1 = param.t1
    t2 = param.t2
    t3 = param.t3
    t4 = param.t4
    t5 = param.t5
    t6 = param.t6
    t7 = param.t7
    t8 = param.t8
    Δxy = param.Δxy
    μ   = param.μ
    ϵ_xz_xz = (kx ::Real, ky ::Real) -> begin
        2*t1*cos(kx) + 2*t2*cos(ky) + 4*t3*cos(kx)*cos(ky) - μ
    end
    ϵ_yz_yz = (kx ::Real, ky ::Real) -> begin
        2*t2*cos(kx) + 2*t1*cos(ky) + 4*t3*cos(kx)*cos(ky) - μ
    end
    ϵ_xy_xy = (kx ::Real, ky ::Real) -> begin
        2*t5*(cos(kx)+cos(ky)) + 4*t6*cos(kx)*cos(ky) + Δxy - μ
    end
    ϵ_xz_yz = (kx ::Real, ky ::Real) -> begin
        4*t4*sin(kx)*sin(ky)
    end
    ϵ_xz_xy = (kx ::Real, ky ::Real) -> begin
        2im*t7*sin(kx) + 4im*t8*sin(kx)*cos(ky)
    end
    ϵ_yz_xy = (kx ::Real, ky ::Real) -> begin
        2im*t7*sin(ky) + 4im*t8*sin(ky)*cos(kx)
    end

    (kx ::Real, ky ::Real) -> begin
        ϵ11 = ϵ_xz_xz(kx,ky)
        ϵ12 = ϵ_xz_yz(kx,ky)
        ϵ13 = ϵ_xz_xy(kx,ky)
        ϵ21 = conj(ϵ12)
        ϵ22 = ϵ_yz_yz(kx,ky)
        ϵ23 = ϵ_yz_xy(kx,ky)
        ϵ31 = conj(ϵ13)
        ϵ32 = conj(ϵ23)
        ϵ33 = ϵ_xy_xy(kx,ky)
        return [ϵ11 ϵ12 ϵ13;
                ϵ21 ϵ22 ϵ23;
                ϵ31 ϵ32 ϵ33]
    end
end

function makeRealspaceDifference(param ::DaghoferParameter, nx ::Integer, ny ::Integer)
    hoppings = hoppingList(param)
    orbitalIndices = Dict("xz"=>1, "yz"=>2, "xy"=>3)
    hamiltonian = zeros(Complex128, nx, ny, 3, 3)
    for (dx, dy, orb1, orb2, val) in hoppings
        iorb1 = orbitalIndices[orb1]
        iorb2 = orbitalIndices[orb2]
        hamiltonian[mod(dx, nx)+1, mod(dy, ny)+1, iorb1, iorb2] += val
        if iorb1 != iorb2
            hamiltonian[mod(-dx, nx)+1, mod(-dy, ny)+1, iorb2, iorb1] += conj(val)
        end
    end
    hamiltonian
end


function makeRealspaceDense(param::DaghoferParameter, nx ::Integer, ny ::Integer, periodic ::Bool = true)
    hoppings = hoppingList(param)
    orbitalIndices = Dict("xz"=>1, "yz"=>2, "xy"=>3)
    hamiltonian = zeros(Complex128, 3, nx, ny, 3, nx, ny)
    for (dx, dy, orb1, orb2, val) in hoppings
        iorb1 = orbitalIndices[orb1]
        iorb2 = orbitalIndices[orb2]
        for ix=0:(nx-1), iy=0:(ny-1)
            jx = ix + dx
            jy = iy + dy
            if periodic
                hamiltonian[iorb1, ix+1, iy+1, iorb2, mod(jx, nx)+1, mod(jy, ny)+1] += val
                if iorb1 != iorb2
                    hamiltonian[iorb2, mod(jx, nx)+1, mod(jy, ny)+1, iorb1, ix+1, iy+1] += conj(val)
                end
            else
                if 0 <= ix < nx && 0 <= iy < ny && 0 <= jx < nx && 0 <= jy < ny
                    hamiltonian[iorb1, ix+1, iy+1, iorb2, jx+1, jy+1] += val
                    if iorb1 != iorb2
                        hamiltonian[iorb2, jx+1, jy+1, iorb1, ix+1, iy+1] += conj(val)
                    end
                end
            end
        end
    end
    hamiltonian
end


function makeRealspaceSparse(param::DaghoferParameter, nx ::Integer, ny ::Integer)
    hoppings = hoppingList(param)
    shape = (3, nx, ny)
    orbitalIndices = Dict("xz"=>1, "yz"=>2, "xy"=>3)
    rows, cols, vals = Int64[], Int64[], Complex128[]

    for (dx, dy, orb1, orb2, val) in hoppings
        iorb1 = orbitalIndices[orb1]
        iorb2 = orbitalIndices[orb2]
        for ix=0:(nx-1), iy=0:(ny-1)
            append!(rows, sub2ind(shape, iorb1, ix+1, iy+1))
            append!(cols, sub2ind(shape, iorb2, mod(ix+dx, nx)+1, mod(iy+dy, ny)+1))
            append!(vals, val)
            if iorb1 != iorb2
                append!(cols, sub2ind(shape, iorb1, ix+1, iy+1))
                append!(rows, sub2ind(shape, iorb2, mod(ix+dx, nx)+1, mod(iy+dy, ny)+1))
                append!(vals, conj(val))
            end
        end
    end
    Hermitian(sparse(rows, cols, vals, 3*nx*ny, 3*nx*ny))
end
=#

end
