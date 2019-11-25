module ToolBox

# Export submodule
export SlcUtil, Load, Geometry, Misc

# export core functions
export SlcRaw



# include submodules from stc folder
include("src/SlcUtil.jl")
include("src/Load.jl")
include("src/Geometry.jl")
include("src/Misc.jl")

using .SlcUtil



end
