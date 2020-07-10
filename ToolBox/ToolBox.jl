module ToolBox

# Export submodule
export SlcUtil, Load, Geometry, Misc, UnitTest

# export core functions
export SlcRaw



# include submodules from stc folder
include("src/SlcUtil.jl")
include("src/Load.jl")
include("src/Geometry.jl")
include("src/Misc.jl")
include("src/UnitTest.jl")

# Why?
using .SlcUtil



end
