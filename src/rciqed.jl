module rciqed

librciqed_so = grasp_path_file = joinpath(@__DIR__, "../deps/libgrasp-rci-qed.so")

function parameter_def(name::Symbol)
    value = Ref{Cint}()
    rv = ccall(
        (:libgrasp_rciqed_parameter_def, librciqed_so), Cint,
        (Cstring, Ref{Cint}),
        name, value,
    )
    rv == 0 || error("Error in libgrasp_rciqed_parameter_def(\"$name\") -> $rv")
    return value[]
end

parameter_def() = NamedTuple(
    name => parameter_def(name)
    for name in (:KEYORB, :NNNP, :NNN1, :NNNW, :NNNWM1, :NNNWM2)
)

end
