

#Peak result type
type JLabPeakset
    summits::Array{Int64,1}
    lbnds::Array{Int64,1}
    rbnds::Array{Int64,1}
    
    JLabPeakset() = new(Array(Int64,0),Array(Int64,0),Array(Int64,0))
end


# #############################################################################
#Naive peak finding methods

#Find peak locations and bounds. Include terraces in bounds
function getmountains{T}(y::Array{T,1};thresh::Union(Real,Array{Real,1})=0)
    
    #Get crude variation
    dy::Array{T,1} = [0,y[2:end]-y[1:end-1]]
    small = (abs(dy) .<= thresh)
    dy[small] = 0; dy[!small] = sign(dy[!small])
    
    #Initiate
    pks = JLabPeakset()
        
    #Initiate
    push!(pks.rbnds,0) #add dummy value
    ascending::Bool = false #start descending
    for i=1:length(y)-1
        if ascending
            #Wait for next drop
            if dy[i+1] == -1
                #Start of descend
                push!(pks.rbnds,i+1)
                ascending = false
            elseif y[i+1] > peakval
                #Ascend continues
                pks.summits[end] = i+1
                peakval = y[i+1]
            end
        else
            #Wait for next rise
            if dy[i+1] == 1
                #Start of new peak
                push!(pks.lbnds,i)
                push!(pks.summits,i+1)
                peakval = y[i+1]
                ascending = true
            elseif dy[i+1] == -1
                #Descend continues
                pks.rbnds[end] = i+1
            end
        end
    end
    
    if ascending
        #Add right bound of last peak
        push!(pks.rbnds,length(y))
    end
    
    #Remove dummy value
    shift!(pks.rbnds)
    
    return pks
end


#Find peak locations and bounds. Exclude terraces in bounds
function gethills{T}(y::Array{T,1};thresh::Union(Real,Array{Real,1})=0)
    
    #Get crude variation
    dy::Array{T,1} = [0,y[2:end]-y[1:end-1]]
    small = (abs(dy) .<= thresh)
    dy[small] = 0; dy[!small] = sign(dy[!small])
    
    #Initiate
    pks = JLabPeakset()
    
    #Initiate
    ascending::Bool = false #start descending
    peakval::T = 0
    for i=1:length(y)-1
        if ascending
            if dy[i+1] == 1 && dy[i] == 0
                #Correct starting foot of new hill
                pks.lbnds[end] = i
            end
            
            #Wait for next drop
            if y[i+1] > peakval
                #Ascend continues
                pks.summits[end] = i+1
                peakval = y[i+1]
            elseif dy[i+1] == -1
                #Start of descend
                ascending = false
            end
        else
            if dy[i+1] >= 0 && length(pks.rbnds) < length(pks.lbnds)
                #Ending foot of hill
                push!(pks.rbnds,i)
            end
            
            #Wait for next rise
            if dy[i+1] == 1
                #Starting foot of new hill
                push!(pks.lbnds,i)
                push!(pks.summits,i+1)
                peakval = y[i+1]
                ascending = true
            end
        end
    end
    
    if length(pks.rbnds) < length(pks.lbnds)
        #Add right bound of last peak
        push!(pks.rbnds,length(y))
    end
    
    return pks
end


# #############################################################################
#Peak processing

#Filter peaks
function filterpeaks{T}(pks::JLabPeakset,y::Array{T,1},ys::Array{Float64,1})
    pks2 = JLabPeakset()
    for i=1:length(pks.summits)
        #Get smoothed height
        altitude = ys[pks.summits[i]]
        leftht = altitude-ys[pks.lbnds[i]]
        rightht = altitude-ys[pks.rbnds[i]]
        ht = 0.5*(leftht+rightht)
        baseht = 0.5*(ys[pks.lbnds[i]]+ys[pks.rbnds[i]])
        
        #Collect isolated pks.summits
        thresh1 = sqrt(maximum([baseht,10]))
        if (maximum([leftht,rightht]) > 4*thresh1) && (minimum([leftht,rightht]) > 2*thresh1)
            #Get unsmoothed height
            altitude = y[pks.summits[i]]
            leftht = altitude-y[pks.lbnds[i]]
            rightht = altitude-y[pks.rbnds[i]]
            ht = ((pks.summits[i]-pks.lbnds[i])*rightht + (pks.rbnds[i]-pks.summits[i])*leftht)/(pks.rbnds[i]-pks.lbnds[i])
            
            if ht > 0
                push!(pks2.summits,pks.summits[i])
                push!(pks2.lbnds,pks.lbnds[i])
                push!(pks2.rbnds,pks.rbnds[i])
            end
        end
    end
    
    return pks2
end

#Fit function to peaks
function fitpeaks{T}(pks::JLabPeakset,x::Array{T,1},y::Array{T,1},ys::Array{Float64,1};peakshape::String="gaussian")
    
    #Initiate
    pkgrps = Array(Array{Int64,1},0)
    fitdata = Array(Array{Float64,1},0)
    fiterrdata = Array(Array{Float64,1},0)
    
    #Fitting functions
    #Individual peak functions take amplitude,width,location as parameters
    if peakshape == "gaussian"
        f_single(x,p) = p[1]*exp(-0.5*(x-p[3]).^2/p[2]^2)
    elseif peakshape == "lorentzian"
        f_single(x,p) = p[1]./(1+(x-p[3]).^2/p[2]^2)
    end
    #Peak group function takes individual peak parameters, then optionally the baseline
    function f_total(x,p)
        if length(p)%3 == 2
            ans = p[end-1]+p[end]*x
        elseif length(p)%3 == 1
            ans = p[end]
        else
            ans = 0
        end
        for k=1:div(length(p),3)
            ans += f_single(x,p[3*k-2:3*k])
        end
        return ans
    end
    
    #Begin
    i = 1
    while i <= length(pks.summits)
        #Collect nearby peaks
        pkgrp = [i]
        while i+1 <= length(pks.summits) && pks.lbnds[i+1]-pks.rbnds[i] < 3
            #next peak is near enough
            i += 1
            push!(pkgrp,i)
        end
        
        #Setup initial parameters
        p0 = Array(Float64,0)
        for j=1:length(pkgrp)
            pkalt = y[pks.summits[pkgrp[j]]]
            leftalt = y[pks.lbnds[pkgrp[j]]]
            rightalt = y[pks.rbnds[pkgrp[j]]]
            pkht = pkalt - 0.5*(leftalt+rightalt)
            push!(p0,pkht) #Estimate peak height
            push!(p0,0.09*(x[pks.rbnds[pkgrp[j]]]-x[pks.lbnds[pkgrp[j]]])) #Estimate peak width
            push!(p0,x[pks.summits[pkgrp[j]]]) #Estimate peak location
        end
        
        push!(pkgrps,pkgrp)
        grprng = [pks.lbnds[pkgrp[1]]:pks.rbnds[pkgrp[end]]]
        
        #Fitting function for multiple peaks
        dof = length(grprng) - (2+3*length(pkgrp))
        if dof > 0
            # grad = (y[grprng[end]]-y[grprng[1]])/(x[grprng[end]]-x[grprng[1]])
            # push!(p0,y[grprng[1]] - grad*x[grprng[1]])
            # push!(p0,grad)
            push!(p0,minimum([y[grprng[1]],y[grprng[end]]]))
            push!(p0,0.0)
        elseif dof == 0
            push!(p0,0.5*(y[grprng[1]]+y[grprng[end]]))
        elseif dof <= -1
            #Reject peak
            push!(fitdata,nans(2*(3*length(pkgrp)+2)))
            i += 1
            continue
        end
        
        #Show parameter guess
        xs = linspace(x[grprng[1]],x[grprng[end]],200)
        #plt.plot(xs,f_total(xs,p0),"k")
        
        #Conduct fit
        nnz = find(ys[grprng] .> 0)
        #yerr = ones(Float64,length(nnz))
        yerr = sqrt(ys[grprng[nnz]])
        fit = (
            try
                getfit(f_total,float64(x[grprng[nnz]]),float64(y[grprng[nnz]]),yerr,p0)
            catch e
                #Reject peak
                show(e)
                println("rejected")
                push!(fitdata,nans(2*(3*length(pkgrp)+2)))
                i += 1
                continue
            end)
        plt.plot(xs,f_total(xs,fit.params),"r")
        
        #Record result
        push!(fitdata,[fit.params,zeros(Float64,2-(length(fit.params)%3))])
        push!(fiterrdata,[fit.perrors,zeros(Float64,2-(length(fit.perrors)%3))])
        
        i += 1
    end
    #plt.show()
    
    return pkgrps,fitdata,fiterrdata
end