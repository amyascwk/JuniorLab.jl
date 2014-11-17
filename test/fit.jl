using Base.Test
using JLab
import PyPlot.plt

# #############################################################################
#Testing "getfit"
println("########## Testing \"getfit\" ##########")

#Sample data
x = linspace(0,10,100)
y = 1 + 3*x
y0 = y + 0.1*randn(100)
y0err = fill(0.1,length(y))
#y1 = y + (0.01*x+0.1).*randn(100)
y1 = [0.960833039028134,1.407268029451956,1.7041362562849804,1.7439408410837736,2.408365436088803,2.3558604670692187,2.6911766100562473,2.985346937471822,3.4325552973746216,3.6518536870040856,3.8198571085348556,4.249638029398669,4.5191875756320465,4.869643189063813,5.200098855165331,5.6028720774344825,5.8669985679112155,6.13446639319173,6.6589983913598445,6.583599635749586,7.0643398456691475,7.22881263766044,7.766208155795485,7.9279248867527095,8.031876149276476,8.622312209495371,8.91270431877012,9.253120846258275,9.46051424213136,9.864708345572803,10.196800366860534,10.262262734504358,10.703336679844554,11.159367909410706,11.448384037553025,11.419407340631926,11.882087020095623,12.123486147998683,12.589137110943783,12.915894183324713,13.195867464426314,13.381066468422823,13.498128075316496,13.946806258004903,14.140720490641081,14.65645590142797,14.92980075431186,15.422243674966104,15.606744487089607,15.664750441796148,16.035908118733797,16.4953299199874,16.958066212026043,16.994285363697312,17.354783816298067,17.700581921115532,18.0306766566351,18.382659996933242,18.231343132691336,19.063032842217314,19.255377262197005,19.81437282587234,19.956484257580982,19.97657224631401,20.360730031763094,20.65882754494313,21.106841540004787,21.252496222358445,21.766760201336467,22.19490679290051,22.6047978108247,22.529362909631008,22.915135996360842,23.17387598157397,23.43698104498548,23.598722726534746,23.73428978953588,24.241320075552377,24.65076910198924,24.78716650508645,25.05422998129448,25.58384488593317,25.93224573773362,25.843980681757863,26.45408533700855,26.888994042344866,27.25415441864504,27.83797148019114,27.550899599039617,28.03717936818957,28.18093241086242,28.548423792990217,28.591866793393848,29.208822396144086,29.49443477158012,29.793649140093873,29.867114283987977,30.526713526135133,30.7660738194201,30.89530630541283]
y1err = 0.01*x + 0.1

#Model
f(x,p) = p[1] + p[2]*x

#Get fit
fit = getfit(f,x,y,ones(Float64,length(y)),[0.,0.])
fit0 = getfit(f,x,y0,y0err,[0.,0.])
fit1 = getfit(f,x,y1,y1err,[0.,0.])

#Fit fields
println("Testing JLabFitresult fields ...")
@test fit0.xdata == x
@test fit0.ydata == y0
@test fit0.yerrors == y0err
@test fit0.model == f
@test fit0.dof == 98

#Fit accuracy
println("Testing accuracy ...")
@test_approx_eq_eps(fit.params[1],1,1e-10)
@test_approx_eq_eps(fit.params[2],3,1e-10)

p0 = sgolaycoeff(x,1)*y0
@test_approx_eq_eps(fit0.params[1],p0[1],1e-7)
@test_approx_eq_eps(fit0.params[2],p0[2],1e-10)

@test_approx_eq_eps(fit1.params[1],1,fit1.perrors[1])
@test_approx_eq_eps(fit1.params[2],3,fit1.perrors[2])

println("Function \"getfit\" passes tests.")