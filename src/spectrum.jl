import PyPlot.plt

# #############################################################################
#Spectrum datafile I/O

#Define spectrum rawdata type
immutable JLabSpectrum
    date::String
    real_time::Int64
    live_time::Int64
    bins::Array{Int64,1}
    counts::Array{Int64,1}
    rois::Array{Int64,2}
end


#Read ascii .spe file
function readspe(filename::String)
    #Read file
    file::IOStream = open(filename,"r")
    lines::Array{String,1} = readlines(file)
    
    #Find fields
    date_mea_idx::Int64 = findfirst(lines,"\$DATE_MEA:\r\n")
    meas_tim_idx::Int64 = findfirst(lines,"\$MEAS_TIM:\r\n")
    data_idx::Int64 = findfirst(lines,"\$DATA:\r\n")
    roi_idx::Int64 = findfirst(lines,"\$ROI:\r\n")
    
    #Extract data -------------------------------------------------------------
    #Date string
    date_mea::String = chomp(lines[date_mea_idx+1])
    #Measurement times (Live, Real)
    meas_tim::Array{Int64,1} = int64(split(lines[meas_tim_idx+1]," "))
    #Data
    binlim::Array{Int64,1} = int64(split(lines[data_idx+1]," "))
    bins::Array{Int64,1} = [binlim[1]:binlim[2]]
    counts::Array{Int64,1} = int64(lines[(data_idx+2):(data_idx+2+binlim[2]-binlim[1])])
    start = 1
    while counts[start] == 0; start += 1; end
    bins = bins[start:end]
    counts = counts[start:end]
    #ROI
    roinum::Int64 = int64(lines[roi_idx+1])
    rois::Array{Int64,2} = Array(Int64,roinum,2)
    for i = 1:roinum
        rois[i,:] = int64(split(lines[roi_idx+1+i]))
    end
    
    return JLabSpectrum(date_mea,meas_tim[1],meas_tim[2],bins,counts,rois)
end


#Show bins
function plotspec(spec::JLabSpectrum;logged::Bool=true)
    if logged
        y::Array{Float64,1} = zeros(Float64,length(spec.counts))
        nonzs::Array{Int64,1} = find(spec.counts)
        y[nonzs] = log(spec.counts[nonzs])
        ylabel::String = "log(Counts)"
    else
        y = float64(spec.counts)
        ylabel = "Counts"
    end
    
    for i=1:length(spec.counts)
        color = [0,0,0.75]
        for j=1:size(spec.rois,1)
            if (spec.rois[j,1] <= spec.bins[i] <= spec.rois[j,2])
                color = [0,0.75,0.75]
                break
            end
        end
        plt.plot([spec.bins[i],spec.bins[i]],[0,y[i]],color=color)
    end
    plt.xlabel("Bin Number")
    plt.ylabel(ylabel)
    plt.axis([spec.bins[1],spec.bins[end],0,1.1*maximum(y)])
    return nothing
end

