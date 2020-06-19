# Implementations of CGD for linear games
using LinearAlgebra
using StaticArrays
using LinearMaps
using IterativeSolvers
using BlockArrays

# The dual of the quadratic distance
# ̄ψ = sum( α .* x.^2 ) / 2
function squared_distance(α) 
    # Computing the primal point associated to the dual point x̄ 
    function Dxψ(x̄)
      return α .\ x̄
    end

    function Dxψ_inv(x̄)
      return α .* x̄
    end

    # Hessian of the dual potential, evaluated in the primal point x
    function Dxxψ(x)
      function Dxxψ_v!(hvp, v)
        hvp .= v ./ α
      end
      return LinearMap{eltype(x)}(Dxxψ_v!, length(x); ismutating=true, issymmetric=true, ishermitian=true) 
    end

    # Hessian of the primal potential, evaluated in the primal point x
    function Dxxψ_inv(x)
      function Dxxψ_v!(hvp, v)
        hvp .= v .* α
      end
      return LinearMap{eltype(x)}(Dxxψ_v!, length(x); ismutating=true, issymmetric=true, ishermitian=true) 
    end
    return Dxψ, Dxψ_inv, Dxxψ, Dxxψ_inv
end

# The dual of the quadratic distance
# ̄ψ = sum( α .* x.^2 ) / 2
function positive_projection(α) 
    # Computing the primal point associated to the dual point x̄ 
    function Dxψ(x̄)
      x̄ .= max.(x̄, zero(eltype(x̄)))
      return α .\ x̄
    end

    function Dxψ_inv(x̄)
      x̄ .= max.(x̄, zero(eltype(x̄)))
      return α .* x̄
    end


    # Hessian of the dual potential, evaluated in the primal point x
    function Dxxψ(x)
      function Dxxψ_v!(hvp, v)
        hvp .= v ./ α
      end
      return LinearMap{eltype(x)}(Dxxψ_v!, length(x); ismutating=true, issymmetric=true, ishermitian=true) 
    end

    # Hessian of the primal potential, evaluated in the primal point x
    function Dxxψ_inv(x)
      function Dxxψ_v!(hvp, v)
        hvp .= v .* α
      end
      return LinearMap{eltype(x)}(Dxxψ_v!, length(x); ismutating=true, issymmetric=true, ishermitian=true) 
    end
    return Dxψ, Dxψ_inv, Dxxψ, Dxxψ_inv
end



# The dual of the (negative) shannon entropy 
# ̄ψ = sum( α .* (x .* log.(x) .- x) )
function shannon_entropy(α) 
    # Computing the primal point associated to the dual point x̄ 
    function Dxψ(x̄)
      return exp.(α .\ x̄)
    end

    function Dxψ_inv(x)
      return α .* log.(x)
    end

    # Hessian of the dual potential, evaluated in the primal point x
    function Dxxψ(x)
      function Dxxψ_v!(hvp, v)
        hvp .= v .* x ./ α
      end
      return LinearMap{eltype(x)}(Dxxψ_v!, length(x); ismutating=true, issymmetric=true, ishermitian=true) 
    end

    # Hessian of the primal potential, evaluated in the primal point x
    function Dxxψ_inv(x)
      function Dxxψ_v!(hvp, v)
        hvp .= v ./ x .* α
      end
      return LinearMap{eltype(x)}(Dxxψ_v!, length(x); ismutating=true, issymmetric=true, ishermitian=true) 
    end
    return Dxψ, Dxψ_inv, Dxxψ, Dxxψ_inv
end

# Combines a list of potentials to act on different parts of a block-vector
function combine_potentials(potentials...)
  function Dxψ(x̄::AbstractBlockVector)
    out = similar(x̄)
    @assert length(potentials) == blocklength(x̄)
    for k = 1 : length(potentials)
      out[Block(k)] .= potentials[k][1](x̄[Block(k)]) 
    end 
    return out
  end

  function Dxψ_inv(x̄::AbstractBlockVector)
    out = similar(x̄)
    @assert length(potentials) == blocklength(x̄)
    for k = 1 : length(potentials)
      out[Block(k)] .= potentials[k][2](x̄[Block(k)]) 
    end 
    return out
  end


  # Hessian of the dual potential, evaluated in the primal point x
  function Dxxψ(x)
    @assert length(potentials) == blocklength(x)
    block_diagonal = [potentials[k][3](x[Block(k)]) for k = 1 : length(potentials)]

    return LinearMaps.BlockDiagonalMap(block_diagonal...)
  end

  # Hessian of the primal potential, evaluated in the primal point x
  function Dxxψ_inv(x)
    @assert length(potentials) == blocklength(x)
    block_diagonal = [potentials[k][4](x[Block(k)]) for k = 1 : length(potentials)]

    return LinearMaps.BlockDiagonalMap(block_diagonal...)
  end

  return Dxψ, Dxψ_inv, Dxxψ, Dxxψ_inv
end


