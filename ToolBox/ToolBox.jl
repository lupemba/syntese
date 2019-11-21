module ToolBox

# Export submodule
export Submod1, Submod2



# export core functions
export hello, fuck_you

# include submodules from stc folder
include("src/Submod1.jl")

# Make submodule in file from src
module Submod2
    export fuck_you, shit

    include("src/fuck_you.jl")
    include("src/shit.jl")

end
using .Submod2

## core functions
function hello(name)
    if is_string(name)
        message = "hello "*name* ". It is nice to meet you"
        return println(message)
    else
        return println("input not a string")
    end

end

function is_string(text)
    return typeof(text)==typeof("String")
end




end
