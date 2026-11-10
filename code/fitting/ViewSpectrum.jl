using SpectralFitting, Plots, LaTeXStrings, ArgParse
pyplot()

# Unit conversions
cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

# Default plotting parameters
default(
    # Fonts
    titlefont = (10, "serif"), 
    guidefont = (8, "serif"), 
    legendfont = (6, "serif"), 
    tickfont = (6, "serif"), 
    colorbar_titlefont = (10, "serif"),
    # Grids
    gridalpha=0.5,
    minorgridalpha=0.2,
    # Image
    dpi=1200,
    size=cm2px.((8,6)),
    framestyle=:box
)

# Data loading function
function LoadData(
    spectrum::String; 
    dataRange::Tuple{Real, Real}=(2,10)
    )

    # Loading in the spectrum file
    data = OGIPDataset(spectrum)

    # Regrouping, normalising, dropping bad channels and curtailing
    regroup!(data)
    # normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data, dataRange...)
end

function ParseArguments()
    s = ArgParseSettings()

    #Add arguments with their default values
    @add_arg_table s begin
        "--OBS_ID"
            help="The observation ID"
            default="MISSING"
        "--Instrument"
            help="XMM EPIC Instrument name (pn, mos1, mos2)"
            default="pn"
        "--INDIR"
            help="The input directory"
            default=pwd()
        "--OUTDIR"
            help="The output directory"
            default=pwd()
        "--XLOG"
            arg_type=Bool
            help="Whether the x axis should use a log scale"
            default=true
        "--YLOG"
            arg_type=Bool
            help="Whether the y axis should use a log scale"
            default=true
    end

    return parse_args(ARGS, s)
end

args = ParseArguments()

# Loading the data (requires that the spectrum files are grouped)
println("Loading Data...")
data = LoadData(
    "$(args["INDIR"])/$(args["OBS_ID"])_$(args["Instrument"])_pi_rebinned.fits"
)

# Setting the observation id
data.user_data.observation_id = args["OBS_ID"]

# Plotting
majorticks = collect(2:10)
vline([6.4]; c=:red, ls=:dash, lw=0.7, label=L"Fe K$\alpha$")

if args["XLOG"]
    xaxis!(:log)
end
if args["YLOG"]
    yaxis!(:log)
end

plt = plot!(
    data; ylabel=L"Counts (s$^{-1}$ keV$^{-1}$)", 
    xticks=(majorticks, majorticks),
    c=:black, markercolor=:black, lc=:black, lw=0.5, msw=0.,
    framestyle=:box, xlims=(minimum(majorticks), maximum(majorticks))
)

# Saving the figure
savefig(
    plt, "$(args["OUTDIR"])/$(args["OBS_ID"])_$(args["Instrument"])_spectrum.png"
)
