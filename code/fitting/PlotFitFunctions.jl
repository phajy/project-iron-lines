using Plots, CSV, DataFrames, Measures

# Unit conversions
cm2in(x) = 0.3937008x
cm2px(x) = Int(round(1.1 * cm2in(x) * 100))

# Setting plotting defaults
default(
    # Fonts
    titlefont = (12, "times"), 
    guidefont = (12, "times"), 
    legendfontfamily = "times",
    legendfontsize = 12,
    legendtitlefontfamily = "times",
    legendtitlefontsize = 12,
    tickfont = (10, "times"), 
    colorbar_titlefont = (10, "times"),
    legendtitlefontvalign = :center,
    # Grids
    gridalpha=0.2,
    minorgridalpha=0.05,
    minorgrid=false,
    minorticks=5,
    # Image
    dpi=400,
    size=1.5*cm2px.([8,7]),
    framestyle=:box
)

function ParseQDP(filename)

    # Open the file
    fileString = open(filename) do file
        read(file, String)
    end

    # Removing the header from the top of the file
    fileString = fileString[findfirst("!", fileString)[end]+2:end]

    # Getting the number of columns
    numcols = count(" ", fileString[1:findfirst("\n", fileString)[1]]) + 1

    # Splitting the datasets
    splitter = join([("NO "^numcols)[1:end-1], "\n"])
    strings = split(fileString, splitter)

    # Generating numbered header
    header = join([join(1:numcols, " "), "\n"])

    # Creating an empty array for each dataset
    output = Array{DataFrame}(undef, length(strings))

    # Assembling dataframes
    for (i, string) in enumerate(strings)
        output[i] = DataFrame(CSV.File(IOBuffer(join([header, String(string)]))))
    end

    return output
end

function PlotDataModel(data; xticks=[2.5, 5, 10])
    # Plotting the data and model
    dataplot = scatter(
        data[!, "1"], data[!, "3"], xerr=data[!, "2"], yerr=data[!, "4"], 
        msw=0.2, ms=0.1, lw=0.5, c=:black, ylabel="counts / s / keV", 
        yaxis=:log, xaxis=:log, xticks=xticks, label="Data", 
        xlims=(first(xticks),last(xticks))
    )
    plot!(dataplot, data[!, "1"], data[!, "5"]; c=:red, label="Model")
    plot!(dataplot, xformatter= _-> "")
end

function PlotRatio(data; xticks=[2.5, 5, 10], yticks=[0.9, 1, 1.1])
    # Plotting the ratio
    ratioplot = hline(
        [1]; c=:green, label=nothing, xlabel="Energy (keV)", ylabel="Ratio", 
        xaxis=:log, xticks=(xticks, string.(xticks)), xlims=(first(xticks),last(xticks)), 
        yticks=(yticks, string.(yticks)), ylims=(first(yticks),last(yticks))
    )
    scatter!(
        ratioplot, data[!, "1"], data[!, "3"]; xerr=data[!, "2"], 
        yerr=data[!, "4"], msw=0, ms=1, lw=0.5, c=:black, label=nothing
    )
end