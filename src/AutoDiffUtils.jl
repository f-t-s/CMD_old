import ForwardDiff
using ForwardDiff: derivative
using Zygote: gradient

function ∇_x(f) 
  function ∇_xf(x, y)
    return first(gradient(x -> f(x,y), x))
  end
  return ∇_xf
end

function ∇_y(f) 
  function ∇_yf(x, y)
    return first(gradient(y -> f(x,y), y))
  end
  return ∇_yf
end

function D_yx(f)
  ∇_yf = ∇_y(f)  
  function D_yxf(x, y, v)
    # return derivative(h ->  ∇_yf(x + h .* v, y), zero(eltype(v))) 
    return ForwardDiff.partials.(∇_yf(ForwardDiff.Dual{Nothing}.(x, v), ForwardDiff.Dual{Nothing}.(y)), 1)
  end
  return D_yxf
end

function D_xy(f)
  ∇_xf = ∇_x(f)  
  function D_xyf(x, y, v)
    # return derivative(h ->  ∇_xf(x, y + h .* v), zero(eltype(v))) 
    return ForwardDiff.partials.(∇_xf(ForwardDiff.Dual{Nothing}.(x), ForwardDiff.Dual{Nothing}.(y, v)), 1) 
  end
  return D_xyf
end