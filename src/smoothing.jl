import LsqFit

# #############################################################################
#Moving average filtering

#Apply filter
function movavgfilt{T}(y::Array{T,1},N::Int64)
    #Usage:
    #   ys = movavgfilt(y,N)
    #       y       Signal to be smoothed
    #       N       Size of window
    #       ys      Smoothed output
    
    #Initiate
    ys::Array{Float64,1} = Array(Float64,length(y))
    
    #Compute mean
    for i=1:length(y)
        s::Int64 = minimum([floor((N-1)/2),i-1,length(y)-i])
        ys[i] = sum(y[i-s:i+s])/(2*s+1)
    end
    
    return ys
end

# #############################################################################
#Savitzsky-Golay filtering

#Get coefficients
function sgolaycoeff{T<:Real}(z::Array{T,1},deg::Int64,diff::Int64=0)
    #Usage: 
    #   C = sgolaycoeff(z,deg)
    #   C = sgolaycoeff(z,deg,diff)
    #       z       Fitting window (eg: [-5:5])
    #       deg     Degree of polynomial for fit
    #       diff    Number of differentiations (must be <= deg, default 0 for 
    #               smoothing)
    #       C       Matrix for transforming window values to polynomial coefficients,
    #               where (i,j)th entry is the multiplier to the jth value for the 
    #               (i-1) order term in the fitting polynomial.
    
    #Check input
    assert(diff <= deg)
    
    #Solve normal equations
    J::Array{T,2} = broadcast(^,z,transpose(0:deg))
    Jt::Array{T,2} = transpose(J)
    C::Array{Float64,2} = inv(Jt*J)*Jt
    
    #Get coefficients for derivatives
    for i = 1:diff
        C = broadcast(*,[1:size(C,1)-1],C[2:end,:])
    end
    
    return C
end

#Apply filter
function sgolayfilt{T<:Real}(y::Array{T,1},N::Int64;deg::Int64=2,diff::Int64=0)
    #Usage:
    #   ys = sgolayfilt(y,N)
    #   ys = sgolayfilt(x,y,N)
    #   ys = sgolayfilt(..., deg=deg)
    #   ys = sgolayfilt(..., diff=diff)
    #       y       Signal to be smoothed
    #       N       Size of window
    #       x       Distribution of signal events (default [1:length(y)] for uniform)
    #       deg     Degree of polynomial for fit
    #       diff    Number of differentiations (must be <= deg, default 0 for 
    #               smoothing)
    #       ys      Smoothed signal
    
    #Check input
    if diff > deg
        error("Differential order cannot exceed polynomial fit degree.")
    end
    m::Int64 = div(N-1,2)
    if 2*m+1 <= deg
        error("Fitting window N = $N is too small for polynomial fit of degree $deg.")
    end
    
    #Initiate
    C::Array{Float64,2} = sgolaycoeff([-m:m],deg,diff)
    ys::Array{Float64,1} = Array(Float64,length(y))
    
    #Starting
    a::Array{Float64,1} = C*y[1:2m+1]
    for i=1:m
        ys[i] = sum(a.*(i-m-1).^[0:length(a)-1])
    end
    
    #Middle
    C1::Array{Float64,1} = C[1,:][:]
    for i=m+1:length(y)-m
        ys[i] = dot(C1,y[i-m:i+m])
    end
    
    #Ending
    a = C*y[end-2m:end]
    for i=length(y)-m+1:length(y)
        ys[i] = sum(a.*(i-length(y)+m).^[0:length(a)-1])
    end
    
    return ys
end

function sgolayfilt{T<:Real}(x::Array{T,1},y::Array{T,1},N::Int64;deg::Int64=2,diff::Int64=0)
    #Check input
    if diff > deg
        error("Differential order cannot exceed polynomial fit degree.")
    end
    m::Int64 = div(N-1,2)
    if 2*m+1 <= deg
        error("Fitting window N = $N is too small for polynomial fit of degree $deg.")
    end
    
    #Initiate
    ys::Array{Float64,1} = Array(Float64,length(y))
    
    #Starting
    C::Array{Float64,2} = sgolaycoeff(x[1:2m+1],deg,diff)
    a::Array{Float64,1} = C*y[1:2m+1]
    for i=1:m
        ys[i] = sum(a.*x[i].^[0:length(a)-1])
    end
    
    #Middle
    for i=m+1:length(y)-m
        C = sgolaycoeff(x[i-m:i+m],deg,diff)
        a = C*y[i-m:i+m]
        ys[i] = sum(a.*x[i].^[0:length(a)-1])
    end
    
    #Ending
    C = sgolaycoeff(x[end-2*m:end],deg,diff)
    a = C*y[end-2*m:end]
    for i=length(y)-m+1:length(y)
        ys[i] = sum(a.*x[i].^[0:length(a)-1])
    end
    
    return ys
end


# #############################################################################
#Robust Local Weighted Scatter-plot Smooth
function robust_lowess{T}(y::Array{T,1},N::Int64;deg::Int64=1,robust::Bool=true)
    #Initiate
    ys::Array{Float64,1} = Array(Float64,length(y))
    w::Array{Float64,1} = Array(Float64,N)
    s::Int64 = floor((N-1)/2)
    
    #Modelling function
    f(x::Array,p::Array{Float64,1}) = sum(broadcast(^,x,transpose([0:length(p)-1]))*p,2)
    f(x::Int64,p::Array{Float64,1}) = sum(p.*x.^[0:length(p)-1])
    
    #Conduct regression
    for i=1:length(y)
        #Get regression region
        if s > i-1
            range::Array{Int64,1} = [1:N]
        elseif s > length(y)-i
            range = [(length(y)-N+1):length(y)]
        else
            range = [(i-s):(i+s)]
        end
        
        #Compute local weights
        w = (1-(abs(range-i)/maximum(abs(range-i))).^3).^3
        
        #Conduct fit
        if length(range) == 1
            #Fit is always exact
            ys[i] = y[i]
        else
            p0::Array{Float64,1} = zeros(Float64,minimum([deg+1,length(range)-1]))
            fit = LsqFit.curve_fit(f,range-i,y[range],w,p0)
            ys[i] = fit.param[1]
        end
    end
    
    if !robust; return ys; end
    
    for iter=1:5
        #Find outliers
        resid::Array{Float64,1} = ys-y
        madev6::Float64 = 6*median(abs(resid))
        inliers::Array{Bool,1} = (abs(resid) .< madev6)
        
        #Compute Robust weights
        rw::Array{Float64,1} = zeros(Float64,length(y))
        rw[inliers] = (1-(resid/madev6).^2).^2
        
        #Conduct regression
        for i=1:length(y)
            #Get regression region
            if s > i-1
                range::Array{Int64,1} = [1:N]
            elseif s > length(y)-i
                range = [(length(y)-N+1):length(y)]
            else
                range = [(i-s):(i+s)]
            end
            
            #Compute local weights
            w = (1-(abs(range-i)/maximum(abs(range-i))).^3).^3
            
            #Conduct fit
            if length(range) == 1
                #Fit is always exact
                ys[i] = y[i]
            else
                p0::Array{Float64,1} = zeros(Float64,minimum([deg+1,length(range)-1]))
                fit = LsqFit.curve_fit(f,range-i,y[range],w.*rw[range],p0)
                ys[i] = fit.param[1]
            end
        end
    end
    
    return ys
end