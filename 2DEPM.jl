using FFTW
using Random
using Plots

γ_dot_arr = 10.0 .^(range(-3.0,stop=1.0,length=30))
L = 4; num_sample = 50
T_leng = 50000; T_max = 5000 # Should be changed
dt = T_max/T_leng

τ = 1.0; τ_el = 1.0; τ_pl = 1.0; σ_Y = 1.0; μ= 1.0
γ_dot_c = σ_Y/(μ*τ)

function get_integral_term(G_ft, γ_pl_dot)
    convol_term = G_ft .* fft(γ_pl_dot)
    integral_term = ifft(convol_term)
    
    return real(integral_term)
end

G_ft = zeros(L,L)
for x in 1:L
    if (x > L/2) qx = L-x
    else qx = x
    end
    for y in 1:L
        if (y > L/2) qy = L-y
        else qy = y
        end
        
        if ((qx==0) && (qy==0)) G_ft[x,y] = 0
        else G_ft[x,y]=-4*(qx^2)*(qy^2)/(qx^2+qy^2)^2
        end
    end
end

p_pl = 1/(τ_pl)*exp(-1/(τ_pl)) # k = 1
p_el = 1/(τ_el)*exp(-1/(τ_el)) # k = 1

n_data =zeros(num_sample, length(γ_dot_arr))
σ_data =zeros(num_sample, length(γ_dot_arr))
println(γ_dot_c)
println("")

σ_result = zeros(length(γ_dot_arr))
n_result = zeros(length(γ_dot_arr))    

t=0
σ=zeros(L,L); n=zeros(L,L); γ_pl_dot=zeros(L,L);
σ_series = zeros(num_sample, T_leng);
σ_series_result = zeros(T_leng)

for num_idx in 1:num_sample
    println("# = ", num_idx)
    for γ_idx in 1:length(γ_dot_arr)

        τ_p = τ_pl*ones(L,L)
        τ_e = τ_el*ones(L,L) + τ_p

        σ_Y = 1.0
        
        σ=zeros(L,L); n=zeros(L,L); γ_pl_dot=zeros(L,L);

        σ_dot = 0.0

        n_avg = 0.0
        t = 0
        γ_dot = γ_dot_arr[γ_idx]

        for t_idx = 1:T_leng
            # σ_series : time-series of σ (inset of Fig. 1).
            if (γ_idx == 7) σ_series[num_idx, t_idx] += sum(σ)/(L*L) end
            γ_pl_dot = n .* σ *(1/(2*μ*τ))
            σ_dot = μ*γ_dot*ones(L,L) + 2*μ*get_integral_term(G_ft, γ_pl_dot)
            σ = σ + dt*σ_dot

            # local state variable
            for x in 1:L
                for y in 1:L
                    if ((t > τ_p[x,y]) && (t < τ_e[x,y]))
                        if ((σ[x,y] >= σ_Y) && (n[x,y] == 0.0))
                            
                            if (rand(Float64) <= p_pl) n[x,y] = 1; n_avg += 1 end 
                        end
                        τ_p[x,y] = (t + τ_pl)
                    end
                    if (t > τ_e[x,y])
                        if ((n[x,y] == 1) && (rand(Float64) <= p_el)) n[x,y] = 0 end
                        τ_e[x,y] = (t + τ_el)
                    end
                end
            end
            t += dt
        end
        σ_data[num_idx, γ_idx] = sum(σ)/(L*L)
        n_data[num_idx, γ_idx] = n_avg/(T_leng*L*L)

    end
end

for γ_idx in 1:length(γ_dot_arr)
    for n_idx in 1:num_sample
        σ_result[γ_idx] += σ_data[n_idx, γ_idx]
        n_result[γ_idx] += n_data[n_idx, γ_idx]
    end
end

for t_idx in 1:T_leng
    for n_idx in 1:num_sample
        σ_series_result[t_idx] += σ_series[n_idx, t_idx]
    end
end

σ_result /= num_sample
n_result /= num_sample
σ_series_result /= num_sample