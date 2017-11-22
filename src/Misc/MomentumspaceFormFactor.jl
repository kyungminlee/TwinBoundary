module MomentumspaceFormFactor

formfactor = Dict(
  "dx2-y2" => (amplitude::Float64=1.0) -> begin
              function (kx, ky)
                  return 2*amplitude*(cos(kx) - cos(ky))
              end
          end,
  "dxy" => (amplitude::Float64=1.0) -> begin
              function (kx, ky)
                  return 4*amplitude*(sin(kx) .* sin(ky))
              end
          end,
  "s±" => (amplitude::Float64=1.0) -> begin
              function (kx, ky)
                  return 4*amplitude*(cos(kx) .* cos(ky))
              end
          end,
  "s*" => (amplitude::Float64=1.0) -> begin
              function (kx, ky)
                  return 2*amplitude*(cos(kx) + cos(ky))
              end
          end,
  "s" => (amplitude::Float64=1.0) -> begin
              function (kx, ky)
                  return amplitude*ones(Float64, size(kx))
              end
          end
)

end
