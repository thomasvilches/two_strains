b = 45:15:90


aux2 = [i for i in b]
dpr = 0.0
for sd in aux2
    run_param_fix([5;10;20;30],1.0,1.0,true,0.95,false,dpr,0.52,0.0,21,sd)
    run_param_fix([5;10;20;30],1.0,1.0,true,0.95,false,dpr,0.52,0.5,21,sd)
    #run_param_fix([5;10;20;30],1.0,1.0,true,true,fm,false,dpr,0.52,1.0,sd)

end

for sd in aux2
    run_param_fix([5;10;20;30],1.0,1.0,true,0.94,false,dpr,0.80,0.0,28,sd)
    run_param_fix([5;10;20;30],1.0,1.0,true,0.94,false,dpr,0.80,0.5,28,sd)
    #run_param_fix([5;10;20;30],1.0,1.0,true,true,fm,false,dpr,0.52,1.0,sd)

end

for sd in aux2
    run_param_fix([5;10;20;30],1.0,1.0,true,0.95,false,dpr,0.52,0.0,21,sd,true)
    run_param_fix([5;10;20;30],1.0,1.0,true,0.95,false,dpr,0.52,0.5,21,sd,true)
    #run_param_fix([5;10;20;30],1.0,1.0,true,true,fm,false,dpr,0.52,1.0,sd)
end

for sd in aux2
    run_param_fix([5;10;20;30],1.0,1.0,true,0.94,false,dpr,0.80,0.0,29,sd,true)
    run_param_fix([5;10;20;30],1.0,1.0,true,0.94,false,dpr,0.80,0.5,28,sd,true)
    #run_param_fix([5;10;20;30],1.0,1.0,true,true,fm,false,dpr,0.52,1.0,sd)

end
