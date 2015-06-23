function cvec(A)
    Ar = real(A); Ai = imag(A)
    return vcat(Ar, Ai)
end

function cmat(A)
    Ar = real(A); Ai = imag(A)
    return hvcat((2,2), Ar, -Ai, Ai, Ar)
end
