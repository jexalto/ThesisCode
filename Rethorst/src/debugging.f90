program debugging
    implicit none
    
    if(isnan(int_Iv))then
        print*, 'Variable:  int_Iv, Integration step: ', m
        print*, 'Lambda: ', lambda
    elseif(isnan(Iv_mulambda))then
        print*, 'Iv_mulambda'
    elseif(isnan(Kv_p))then
        print*, 'Kv_p'
    elseif(isnan(Iv))then
        print*, 'Iv'
    elseif(isnan(Kv))then
        print*, 'Kv'
    endif

    if(isnan(multiplier))then
        print*, 'multiplier ', Iv_mulambda, Kv_p, Kv, lambda, mu, mulambda
    endif
end program debugging