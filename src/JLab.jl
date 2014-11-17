module JLab
    
    export to_sf, to_dp
    export JLabFitresult, getfit, plotfit
    export movavgfilt, sgolaycoeff, sgolayfilt
    export JLabPeakset, getmountains, gethills, filterpeaks, fitpeaks
    export JLabSpectrum, readspe, plotspec
    
    include("sfdp.jl")
    include("fit.jl")
    include("smoothing.jl")
    include("peaks.jl")
    include("spectrum.jl")
    
end
