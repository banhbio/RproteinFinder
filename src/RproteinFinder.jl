module RproteinFinder
using BlastLCA
using Taxonomy
using SQLite
using FASTX
using Statistics

include("data.jl")
include("externalcommand.jl")
include("build.jl")
include("find.jl")

end 
