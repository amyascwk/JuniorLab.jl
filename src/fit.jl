import LsqFit
import Distributions
import PyPlot.plt

#Fit result type
immutable JLabFitresult
    #Fit inputs
    xdata::Array{Float64,1}     #independent variable
    ydata::Array{Float64,1}     #dependent variable
    yerrors::Array{Float64,1}   #uncertainties in dependent variable
    model::Function             #model to be fitted against
    
    #Fit results
    residue::Array{Float64,1}   #weighted residue
    params::Array{Float64,1}    #best fit parameters for model
    perrors::Array{Float64,1}   #uncertainties in best fit parameters
    dof::Int64                  #number of degrees of freedom in fit
    chisq::Float64              #chi^2 value of fit
    probability::Float64        #probability of chi^2 value
end

#Fitting procedure
function getfit{T1<:Real,T2<:Real,T3<:Real}(model::Function,xdata::Array{T1,1},ydata::Array{T2,1},initparams::Array{T3,1};kwargs...)
    return getfit(model,xdata,ydata,ones(Float64,length(ydata)),initparams;kwargs...)
end

function getfit{T1<:Real,T2<:Real,T3<:Real,T4<:Real}(model::Function,xdata::Array{T1,1},ydata::Array{T2,1},yerrors::Array{T3,1},initparams::Array{T4,1};alpha=0.6827)
    #Takes a model function with experimental data and returns a fit result.
    #model is a function f(x,p) with x an array of values for the independent 
    #variable, and p an array of parameters for the model, and returns values for the 
    #dependent variable.
    #
    #Usage:
    #   fit = getfit(model,xdata,ydata,yerrors,initparams)
    #       model           Function on (x,params) 
    #       xdata           independent variable values
    #       ydata           dependent variable values
    #       yerrors         uncertainty in dependent variable
    #       initparams      initial guess for model parameters
    
    #Conduct fit
    weights = 1./yerrors
    fit = LsqFit.curve_fit(model,float64(xdata),float64(ydata),float64(weights),float64(initparams))
    perrors = LsqFit.estimate_errors(fit,alpha)
    
    chisq = sum(fit.resid.^2)
    chisqdist = Distributions.Chisq(fit.dof)
    probability = Distributions.ccdf(chisqdist,chisq)
    
    return JLabFitresult(float64(xdata),float64(ydata),float64(yerrors),model,fit.resid,fit.param,perrors,fit.dof,chisq,probability)
end

#Plot fit result
function plotfit(fit::JLabFitresult;datakwargs=(),modelkwargs=())
    #Makes a plot of data points and fit curve with annotations on fit parameters
    
    #Plot experimental data
    datapts = plt.errorbar(fit.xdata,fit.ydata,fit.yerrors,fmt=".",label="Expt data";datakwargs...)
    
    #Plot Model
    x = linspace(minimum(fit.xdata),maximum(fit.xdata),200)
    y = fit.model(x,fit.params)
    modelline = plt.plot(x,y,"r",lw=2,label="Model";modelkwargs...)
    
    #Annotations
    println("@sprintf(\"\\\$\\\\chi^2\\\$ / ndf = %.2f / %d = %.2f\\np(\\\$\\\\chi^2\\\$) = %.3f\\n\",fit.chisq,fit.dof,fit.chisq/fit.dof,fit.probability)")
    @printf("\$\\chi^2\$ / ndf = %.2f / %d = %.2f\np(\$\\chi^2\$) = %.3f\n\n",fit.chisq,fit.dof,fit.chisq/fit.dof,fit.probability)
    annot = plt.text(0.5*(minimum(fit.xdata)+maximum(fit.xdata)),0.5*(minimum(fit.ydata)+maximum(fit.ydata)),@sprintf("\$\\chi^2\$ / ndf = %.2f / %d = %.2f\np(\$\\chi^2\$) = %.3f\n\n",fit.chisq,fit.dof,fit.chisq/fit.dof,fit.probability))
    
    #Plot
    plt.xlabel("\$x\$")
    plt.ylabel("\$y\$")
    plt.title("\$y\$ vs \$x\$")
    plt.legend()
    
    return datapts,modelline,annot
end

# #############################################################################
#Common fitting functions

#Linear
linear{T<:Real,N}(x::Array{T,N},p::Array{Float64,1}) = p[1]+p[2]*x

#Polynomial
polynomial{T<:Real,N}(x::Array{T,N},p::Array{Float64,1}) = broadcast(^,x,transpose([0:length(p)-1]))*p

#Gaussian
gaussian{T<:Real,N}(x::Array{T,N},p::Array{Float64,1}=[1/sqrt(2*pi),1.,0.]) = p[1]*exp(-0.5*(x-p[3]).^2/p[2]^2)

#Lorentzian
lorentzian{T<:Real,N}(x::Array{T,N},p::Array{Float64,1}=[2/pi,1.,0.]) = p[1]/(1+(x-p[3]).^2/p[2]^2)

#Exponential decay
decay{T<:Real,N}(x::Array{T,N},p::Array{Float64,1}=[1.,1.]) = p[1]*exp(-p[2]*x)