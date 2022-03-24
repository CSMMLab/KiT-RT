using DelimitedFiles
using PyPlot

# command is julia plot.jl inputfilename outputfilename

println(ARGS)

if length(ARGS) > 1
    input1 = ARGS[1];
    input2 = ARGS[2];
else
    println("ERROR: need to specify input and output file in command line")
end

function Vec2Mat(nx,ny,v)
    m = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            m[i,j] = v[(i-1)*ny + j]
        end
    end
    return m;
end


result1 = readdlm("../../result/$(input1).txt", Float64)
x = result1[:,1];
y = result1[:,2];
nx1 = Int(floor(sqrt(length(x))))
xGrid1 = y[1:nx1]
dose1 = Vec2Mat(nx1,nx1,result1[:,3]);

X1 = (xGrid1'.*ones(nx1))
Y1 = (xGrid1'.*ones(nx1))'

result2 = readdlm("../../result/$(input2).txt", Float64)
x = result2[:,1];
y = result2[:,2];
nx2 = Int(floor(sqrt(length(x))))
xGrid2 = y[1:nx2]
dose2 = Vec2Mat(nx2,nx2,result2[:,3]);

println("maximums are: left ",maximum(dose1),"; right ",maximum(dose2))

X2 = (xGrid2'.*ones(nx2))
Y2 = (xGrid2'.*ones(nx2))'

fig, (ax1, ax2,ax3) = plt.subplots(1, 3,figsize=(30,10),dpi=100)
ax1.pcolormesh(Y1,X1,dose1,vmin=0.0,cmap="magma")
ax2.pcolormesh(Y2,X2,dose2,vmin=0.0,cmap="magma")
ax3.plot(xGrid1,dose1[:,Int(floor(length(xGrid1)/2))]./maximum(dose1[:,Int(floor(length(xGrid1)/2))]), "b--", linewidth=2, label="left, $(input1)", alpha=0.8)
ax3.plot(xGrid2,dose2[:,Int(floor(length(xGrid2)/2))]./maximum(dose2[:,Int(floor(length(xGrid2)/2))]), "r-", linewidth=2, label="right, $(input2)", alpha=0.8)
ax3.legend(loc="upper left")
ax3.set_xlim([xGrid1[1],xGrid1[end]])
#ax.set_ylim([0,1.05])
ax1.tick_params("both",labelsize=20) 
ax2.tick_params("both",labelsize=20) 
ax3.tick_params("both",labelsize=20) 
tight_layout()
savefig("comparison.png")
