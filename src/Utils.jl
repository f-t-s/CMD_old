using BlockArrays

function _concatenate_functions(funcs...)
    function output_function(v)
        out_v = similar(v)
        for (func_ind, func) in enumerate(funcs)
            out_v[Block(func_ind)] = func(v[Block(func_ind)])
        end
        return out_v
    end
end

