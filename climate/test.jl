using ModelingToolkit

@variables t 

function f end

@register_symbolic f(t)

function test()
    global f(t) = t > 0 ? 1 : 0
    return f
end

test()(t)