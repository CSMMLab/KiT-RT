using DelimitedFiles
using PyPlot

function Vec2Mat(nx,ny,v)
    m = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            m[i,j] = v[(i-1)*ny + j]
        end
    end
    return m;
end


result = readdlm("../../result/CSD_1STOrder_HighDensity.txt", Float64)
x = result[:,1];
y = result[:,2];
nx = Int(floor(sqrt(length(x))))
xGrid = y[1:nx]
dose = Vec2Mat(nx,nx,result[:,3]);

X = (xGrid'.*ones(nx))
Y = (xGrid'.*ones(nx))'

fig = figure("Dose, KiT-RT",figsize=(10,10),dpi=100)
ax = gca()
pcolormesh(Y,X,dose,vmin=0.0)
ax.tick_params("both",labelsize=20) 
#colorbar()
plt.xlabel("x", fontsize=20)
plt.ylabel("y", fontsize=20)
plt.title(L"dose, $S_N$", fontsize=25)
tight_layout()
savefig("dose.png")