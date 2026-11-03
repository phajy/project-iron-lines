include("PlotFit.jl")

# Reading the data
data, ratio = ParseQDP("")

xticks=[3,5,10]
yticks=[0.9,1,1.1,1.15]

# Combining the plots
layout = @layout [a{0.7h}; b{0.3h}]
plot(
    PlotDataModel(data; xticks=xticks), 
    PlotRatio(ratio; xticks=xticks, yticks=yticks); 
    layout = layout, bottom_margin=[-7mm 2mm]
)