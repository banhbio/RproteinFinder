module RproteinFinder
using BlastLCA
using Taxonomy
using FASTX
using JSON

include("data.jl")
include("utils.jl")
include("externalcommand.jl")
include("kofamscan.jl")
include("build.jl")
include("find.jl")

end 
