using LinearAlgebra


# Struct containing the information for cmd 
struct ExtraMD 
    # First player
    # ############################################################
    # The present optimization variable for x
    x̄ 

    # An initial guess to use when solving for x
    x̄_guess 

    # The dual potential's derivative that is used to map a dual variable to its primal
    Dxψ

    # The dual potential's Hessian
    Dxxψ 

    # The inverse of the dual potential's Hessian, which is used precondition the conjugate gradient solver 
    Dxxψ_inv

    # Gradient of x's loss function
    ∇_xf

    # Second player
    # ############################################################
    # The present optimization variable for y
    ȳ 

    # An initial guess to use when solving for y
    ȳ_guess 

    # The dual potential's derivative that is used to map a dual variable to its primal
    Dyψ

    # The dual potential's Hessian
    Dyyψ 

    # The inverse of the dual potential's Hessian, which is used precondition the conjugate gradient solver 
    Dyyψ_inv

    # Gradient of ys loss function
    ∇_yg
end

function ExtraMD(x, xPotential, f, y, yPotential, g)
  return ExtraMD(x, zero(x), xPotential[1], xPotential[3], xPotential[4], ∇_x(f), y, zero(y), yPotential[1], yPotential[3], yPotential[4], ∇_y(g))
end

function step!(extramd::ExtraMD)
  x = extramd.Dxψ(extramd.x̄)
  Dxxψ = extramd.Dxxψ(x)
  Dxxψ_inv = extramd.Dxxψ_inv(x)

  y = extramd.Dyψ(extramd.ȳ)
  Dyyψ = extramd.Dyyψ(y)
  Dyyψ_inv = extramd.Dyyψ_inv(y)

  ∇_xf = extramd.∇_xf(x, y)
  ∇_yg = extramd.∇_yg(x, y)

  extramd.x̄_guess .= - Dxxψ_inv * (Dxxψ * ∇_xf)
  extramd.ȳ_guess .= - Dyyψ_inv * (Dyyψ * ∇_yg)


  # Computing the position of the extra step
  extra_x = extramd.Dxψ(extramd.x̄ + extramd.x̄_guess)
  extra_y = extramd.Dxψ(extramd.ȳ + extramd.ȳ_guess)

  # Computing the gradient in the extra step 
  ∇_xf = extramd.∇_xf(extra_x, extra_y)
  ∇_yg = extramd.∇_yg(extra_x, extra_y)

  extramd.x̄_guess .= - Dxxψ_inv * (Dxxψ * ∇_xf)
  extramd.ȳ_guess .= - Dyyψ_inv * (Dyyψ * ∇_yg)

  extramd.x̄ .+= extramd.x̄_guess
  extramd.ȳ .+= extramd.ȳ_guess

  return 4
end

function iterate_primal(extramd::ExtraMD)
  return extramd.Dxψ(extramd.x̄), extramd.Dyψ(extramd.ȳ)
end

function iterate_dual(extramd::ExtraMD)
  return copy(extramd.x̄), copy(extramd.ȳ)
end