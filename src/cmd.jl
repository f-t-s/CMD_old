# Implementations of CGD for linear games
using LinearAlgebra
using StaticArrays
using LinearMaps
using IterativeSolvers
include("AutoDiffUtils.jl")

function testinv!(x, A, b)
  x .= Matrix(A)\b
end

# Struct containing the information for cmd 
struct CMD 
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

    # Mixed Hessian of x's loss function
    D_xyf 

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

    # Mixed Hessian of ys loss function
    D_yxg

    inv_x::Vector{Bool}

    algorithm!
end

function CMD(x, xPotential, f, y, yPotential, g, algorithm=gmres!)
  return CMD(x, zero(x), xPotential[1], xPotential[3], xPotential[4], ∇_x(f), D_xy(f), y, zero(y), yPotential[1], yPotential[3], yPotential[4], ∇_y(g), D_yx(g), trues(1), algorithm)
end

# preconditioner class 
struct Prec
  inverse
end

import LinearAlgebra.ldiv!
function ldiv!(y, P::Prec, x)
  y .= P.inverse * x
end

function ldiv!(P::Prec, x)
  return P.inverse * x
end

function step!(cmd::CMD)
  x = cmd.Dxψ(cmd.x̄)
  Dxxψ = cmd.Dxxψ(x)
  Dxxψ_inv = cmd.Dxxψ_inv(x)

  y = cmd.Dyψ(cmd.ȳ)
  Dyyψ = cmd.Dyyψ(y)
  Dyyψ_inv = cmd.Dyyψ_inv(y)

  function D_xyf!(u, v)
    u .= cmd.D_xyf(x, y, v)
  end

  function D_yxg!(u, v)
    u .= cmd.D_yxg(x, y, v)
  end


  ∇_xf = cmd.∇_xf(x, y)
  ∇_yg = cmd.∇_yg(x, y)

  # wrapping the mutation hvp into linear maps
  D_xyf = LinearMap{eltype(cmd.ȳ)}(D_xyf!, length(cmd.x̄), length(cmd.ȳ); ismutating=true, issymmetric=false, ishermitian=false) 

  D_yxg = LinearMap{eltype(cmd.x̄)}(D_yxg!, length(cmd.ȳ), length(cmd.x̄); ismutating=true, issymmetric=false, ishermitian=false) 

  # For the exact solve

  # cmd.x̄ .+= - Matrix(Dxxψ - Dxxψ * D_xyf * Dyyψ * D_yxg * Dxxψ) \ (Dxxψ * ∇_xf - Dxxψ * D_xyf * Dyyψ * ∇_yg)
 
  # cmd.ȳ .+= - Matrix(Dyyψ - Dyyψ * D_yxg * Dxxψ * D_xyf * Dyyψ) \ (Dyyψ * ∇_yg - Dyyψ * D_yxg * Dxxψ *∇_xf)

  # For timing comparison
 
  if cmd.inv_x[1] == true
    P = Prec(Dxxψ_inv)
    # P = cholesky(Matrix(Dxxψ))
    hist = cmd.algorithm!(cmd.x̄_guess, Dxxψ - Dxxψ * D_xyf * Dyyψ * D_yxg * Dxxψ, -(Dxxψ * ∇_xf - Dxxψ * D_xyf * Dyyψ * ∇_yg), log=true, 
    tol=1e-9,
    Pl=P,
    )[2]
    cmd.ȳ_guess .= - ∇_yg - D_yxg * Dxxψ * cmd.x̄_guess

    cmd.inv_x[1] = false
  else 
    P = Prec(Dyyψ_inv)
    # P = cholesky(Matrix(Dyyψ))
    hist = cmd.algorithm!(cmd.ȳ_guess, Dyyψ - Dyyψ * D_yxg * Dxxψ * D_xyf * Dyyψ, -(Dyyψ * ∇_yg - Dyyψ * D_yxg * Dxxψ *∇_xf), log=true, 
    tol=1e-9,
    Pl=P,
    )[2]
    cmd.x̄_guess .= - ∇_xf - D_xyf * Dyyψ * cmd.ȳ_guess

    cmd.inv_x[1] = true 
  end


  cmd.x̄ .+= cmd.x̄_guess
  cmd.ȳ .+= cmd.ȳ_guess

  return 4 + 2 * hist.iters
end

function iterate_primal(cmd::CMD)
  return cmd.Dxψ(cmd.x̄), cmd.Dyψ(cmd.ȳ)
end

function iterate_dual(cmd::CMD)
  return copy(cmd.x̄), copy(cmd.ȳ)
end





