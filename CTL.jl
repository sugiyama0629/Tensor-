function ising_transfer_element_2D(i1::Int,i2::Int,j1::Int,j2::Int,
    β::Float64, J::Float64)
    σ1 = 2*i1 -1 #left spin
    σ2 = 2*i2 -1 #right spin
    s1 = 2*j1 -1 #down spin
    s2 = 2*j2 -1 #up spin
    return exp(β*J*(σ1*s1 + s1*σ2 + σ2*s2 + s2*σ1))
end
#2Dイジングの転送行列を定義する関数
function ising_transfer_tensor_2D(β::Float64, J::Float64,
    σ1::Index, σ2::Index, s1::Index, s2::Index)
    T = ITensor(σ1, σ2, s1, s2)
    for i1 in 0:1, i2 in 0:1,
        j1 in 0:1, j2 in 0:1
        T[σ1=>i1+1, σ2=>i2+1, s1=>j1+1, s2=>j2+1] =
            ising_transfer_element_2D(i1, i2, j1, j2, β, J)
    end
    return T
end

function fixend(ind1::Index)
    D = ITensor(ind1)
    D[ind1 => 1] = 1.0
    D[ind1 => 2] = 0.0
    return D
end 