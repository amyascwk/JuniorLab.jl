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
    dof::Int64;                 #number of degrees of freedom in fit
    chisq::Float64              #chi^2 value of fit
    probability::Float64        #probability of chi^2 value
end

#Fitting procedure
function getfit(model::Function,xdata::Array{Float64,1},ydata::Array{Float64,1},yerrors::Array{Float64,1},initparams::Array{Float64,1})
    #Takes a model function with experimental data and returns a fit result.
    #model is a function f(x,p) with x an array of values for the independent 
    #variable, and p an array of parameters for the model, and returns values for the 
    #dependent variable.
    #
    #Usage:
    #   fit = getfit(model,xdata,ydata,yerrors,initparams)
    #       model           Function on (x,params) 
    #       xdata           independent variable
    #       ydata           dependent variable
    #       yerrors         uncertainty in dependent variable
    #       initparams      initial guess for model parameters
    
    #Conduct fit
    weights = 1./yerrors;
    fit = LsqFit.curve_fit(model,xdata,ydata,weights,initparams)
    perrors = LsqFit.estimate_errors(fit);
    
    chisq = sum(fit.resid.^2);
    chisqdist = Distributions.Chisq(fit.dof);
    probability = Distributions.ccdf(chisqdist,chisq)
    
    return JLabFitresult(xdata,ydata,yerrors,model,fit.resid,fit.param,perrors,fit.dof,chisq,probability);
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
    
    return datapts,modelline
end