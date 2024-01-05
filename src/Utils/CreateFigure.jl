using Plots, XLSX, DataFrames, Images, ImageView
using OffsetArrays 

function plot_iteration_scalar(Sv::Vector{Float64}, ylab::String; xlab="Iteration")

    mx =maximum(Sv)*1.2
    mn = minimum(Sv)*0.8
    N = length(Sv)
    
    plt_s = plot([0],[Sv[1]],seriestype=:scatter,markercolor=:red,label=false)
    plot!([0],[Sv[1]],linecolor=:red,linewidth=2,label=false)
    plot!(xlims=([0,N]),ylims=([mn,mx]))
    plot!(ylabel=ylab,xlabel=xlab)
    
    savefig(plt_s, joinpath("Figures","plt-0"))
    for i =2:1:N
        plot!(0:i-1,Sv[1:i],seriestype=:scatter,markercolor=:red,label=false)
        plot!(0:i-1,Sv[1:i],linecolor=:red,linewidth=2,label=false)
        savefig(plt_s, joinpath("Figures","plt-$(i-1)"))
    end
    
end


for i = 0:1:9
img1 = load("Figures/C.000$i.png")
img2 = load("Figures/plt-$i.png")

img1_o   = OffsetArray(img1, 0, 0)
img2_o   = OffsetArray(img2,  100,  900)



r,g = paddedviews(zero(eltype(img1_o)), img1_o, img2_o)

out_over = copy(r)
out_over[axes(img2_o)...] .= img2_o

imres = mosaicview(out_over)
save("Figures/res-$i.png", imres)

end
