module RproteinFinder
using BlastLCA
using Taxonomy
using FASTX

include("data.jl")
include("externalcommand.jl")
include("kofamscan.jl")
include("build.jl")
include("find.jl")

include("cli.jl")

end 
