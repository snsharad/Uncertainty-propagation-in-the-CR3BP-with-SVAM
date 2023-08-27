function scaled = xform_uniform(tmpA, tmpB, x, id)

    if id == 1 % original to [-1,1]

        tmpC = x - tmpA;
        tmpD = tmpB .^ -1;
        scaled = tmpD .* tmpC;  
    
    else
        
        scaled = tmpA + (tmpB .* x); %[-1,1] to original
         
    end
    
end