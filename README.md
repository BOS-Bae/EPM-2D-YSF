# Elastoplastic model (EPM)

This is the Julia simulation code for elastoplastic model(EPM) of two dimensional yield stress fluids(YSF).

The code is for reproducing Fig. 1. in paper regarding EPM of YSF : https://journals.aps.org/pre/abstract/10.1103/PhysRevE.71.010501.

### $L=4$

We can check the code with some figures below :

```julia
using LaTeXStrings
plot(γ_dot_arr, ones(length(γ_dot_arr)),xaxis=:log,yaxis=:log, ylims = (0.1,50), label= L"\sigma=\sigma_Y",ls=:dash,lc=:green, xlabel = L" \dot \gamma /\dot \gamma_c ", ylabel = L" \sigma/\sigma_Y " )
scatter!(γ_dot_arr, σ_result,xaxis=:log,yaxis=:log, ylims = (0.1,50), label= L"L=4",ms=2,mc=:blue, xlabel = L" \dot \gamma /\dot \gamma_c ", ylabel = L" \sigma/\sigma_Y " )
```
<img src="https://github.com/BOS-Bae/EPM-2D-YSF/blob/main/Fig1.png" width="400" height="250"/>

```julia
using LaTeXStrings
plot(dt*(1:T_leng)[1:8000], σ_series_result[1:8000], label= L"L=4",xlabel = L" t/ \tau ", ylabel = L" \sigma/\sigma_Y " ) 
```
<img src="https://github.com/BOS-Bae/EPM-2D-YSF/blob/main/Fig1_inset.png" width="400" height="250"/>