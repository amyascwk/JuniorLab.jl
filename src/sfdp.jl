#Get mantissa and exponent parts of a number
function splitnum(x::Real)
    parts::Array{String,1} = split(string(float64(x)),"e")
    if length(parts) == 1
        if x != 0
            exponent::Int16 = floor(log10(abs(parsefloat(parts[1]))))
            mantissa::Float64 = parsefloat("$(parts[1])e$(-exponent)")
        else
            mantissa = 0.0
            exponent = 0
        end
    else
        mantissa = parsefloat(parts[1])
        exponent = parseint(parts[2])
    end
    return mantissa, exponent
end

#Format to n significant figures
function to_sf(x::Real,num_sf::Int64;expformat::String="e%",force=false)
    #Usage:
    #   s = to_sf(x,nsf)
    #   s = to_sf(..., force=force)
    #   s = to_sf(..., expformat=expformat)
    #       x           real number to be converted
    #       nsf         number of significant figures to be presented in
    #       force       true if scientific notation is to be enforced for small 
    #                   exponents in [-1,3] (default false to allow casual notation 
    #                   sometimes)
    #       expformat   string for exponent to be displayed, occurrence of '%' denotes 
    #                   position of exponent value
    #       s           output string of formatted number
    
    #Check input
    if num_sf < 1
        error("Number of significant figures must be an integer greater or equal to 1.")
    end
    
    #Account for zeros
    if x == 0; return "0"; end
    
    #Split number into mantissa and exponent parts
    mantissa::Float64, exponent::Int16 = splitnum(x)
    #Round to correct number of significant figures
    digits::String = string(int64(round(parsefloat("$(mantissa)e$(num_sf-1)"))))
    
    if !force
        #Allow standard notation for...
        if exponent == -1;
            #... numbers in [0.1,1), and...
            return "0."*digits
        elseif num_sf-1 <= exponent <= 3
            #... integers in [1,9999]
            return rpad(digits,exponent+1,"0")
        end
    end
    
    #Use scientific notation
    if length(digits) == 1
        #no decimal point for single digit
        numstr::String = digits
    else
        #Include decimal point
        numstr = "$(digits[1])."*digits[2:end]
    end
    if exponent != 0
        #append exponent
        expsymb = replace(expformat,"%",string(exponent))
        numstr *= expsymb
    end
    
    return numstr
end

#Format to n decimal places
function to_dp(x::Real,num_dp::Int64)
    #Usage:
    #   s = to_dp(x,ndp)
    #       x       real number to be converted
    #       ndp     number of decimal places to be presented in. Can be negative 
    #               to round to places left of decimal point
    #       s       output string of formatted number
    
    if num_dp > 0
        return eval(parse("@sprintf(\"%.$(num_dp)f\",$(float64(x)))"))
    else
        #Get significant digits
        digits::String = string(int64(round(parsefloat("$(float64(x))e$(num_dp)"))))
        if digits == "0"
            return digits
        else
            #pad with zeros
            return digits*"0"^abs(num_dp)
        end
    end
end
