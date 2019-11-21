module Submod1

export hi


function hi(n)
    if is_int(n)
        message = "hi"^n
        return println(message)
    else
        return println("input not a int")
    end

end


function test3(name)
    return println("For the shit of it "* name)
end


function is_int(number)
    return typeof(number)==typeof(1)
end


end
