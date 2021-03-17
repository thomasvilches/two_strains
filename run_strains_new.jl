#run_param_fix(herd_im_v = [0],fs=0.0,fm=0.0,insert_sec=false,strain2_trans=1.0,vaccinate = false,v_e = 0.0,vfd = v_e/2.0,rd=0.0,sdd=21,extra_d=0,ed=false,nsims=1000)

#run_param_fix([10;20],1.0,1.0)##no vaccine one strain
run_param_fix([10],1.0,1.0,true,1.06)#2 strains 10%
run_param_fix([10],1.0,1.0,true,1.155)#2 strains 20%
run_param_fix([10],1.0,1.0,true,1.25)#2 strains 30%
run_param_fix([10],1.0,1.0,true,1.36)#2 strains 40%
run_param_fix([10],1.0,1.0,true,1.47)#2 strains 50%
run_param_fix([10],1.0,1.0,true,1.67)#2 strains 50%

run_param_fix([20],1.0,1.0,true,1.11)#2 strains 10%
run_param_fix([20],1.0,1.0,true,1.2)#2 strains 20%
run_param_fix([20],1.0,1.0,true,1.327)#2 strains 30%
run_param_fix([20],1.0,1.0,true,1.417)#2 strains 40%
run_param_fix([20],1.0,1.0,true,1.494)#2 strains 50%
run_param_fix([20],1.0,1.0,true,1.66)#2 strains 50%
 
for b in [1.06;1.155;1.25;1.36;1.47;1.67]#second strain transmissibility

    ###Pfizer
    
    #50% reduction
     run_param_fix([10],1.0,1.0,true,b,true,"pfizer",0,false)#no extra dose at day 80
    run_param_fix([10],1.0,1.0,true,b,true,"pfizer",15,false)#15 extra doses, keeping priority groups
    run_param_fix([10],1.0,1.0,true,b,true,"pfizer",15,true)#15 extra doses given to general population
    run_param_fix([10],1.0,1.0,true,b,true,"pfizer",30,false)#15 extra doses, keeping priority groups
    run_param_fix([10],1.0,1.0,true,b,true,"pfizer",30,true)#15 extra doses given to general population
    run_param_fix([10],1.0,1.0,true,b,true,"pfizer",45,false)#15 extra doses, keeping priority groups
    run_param_fix([10],1.0,1.0,true,b,true,"pfizer",45,true)#15 extra doses given to general population
    run_param_fix([10],1.0,1.0,true,b,true,"pfizer",60,false)#15 extra doses, keeping priority groups
    run_param_fix([10],1.0,1.0,true,b,true,"pfizer",60,true)#15 extra doses given to general population

    #50% reduction
     run_param_fix([10],1.0,1.0,true,b,true,"moderna",0,false)#no extra dose at day 80
    run_param_fix([10],1.0,1.0,true,b,true,"moderna",15,false)#15 extra doses, keeping priority groups
    run_param_fix([10],1.0,1.0,true,b,true,"moderna",15,true)#15 extra doses given to general population
    run_param_fix([10],1.0,1.0,true,b,true,"moderna",30,false)#30 extra doses, keeping priority groups
    run_param_fix([10],1.0,1.0,true,b,true,"moderna",30,true)#30 extra doses given to general population
    run_param_fix([10],1.0,1.0,true,b,true,"moderna",45,false)#45 extra doses, keeping priority groups
    run_param_fix([10],1.0,1.0,true,b,true,"moderna",45,true)#45 extra doses given to general population
    run_param_fix([10],1.0,1.0,true,b,true,"moderna",60,false)#60 extra doses, keeping priority groups
    run_param_fix([10],1.0,1.0,true,b,true,"moderna",60,true)#60 extra doses given to general population

end


for b in [1.11;1.2;1.327;1.417;1.494;1.66]#second strain transmissibility

    ###Pfizer
    #no reduction
   
    #50% reduction
     run_param_fix([20],1.0,1.0,true,b,true,"pfizer",0,false)#no extra dose at day 80
    run_param_fix([20],1.0,1.0,true,b,true,"pfizer",15,false)#15 extra doses, keeping priority groups
    run_param_fix([20],1.0,1.0,true,b,true,"pfizer",15,true)#15 extra doses given to general population
    run_param_fix([20],1.0,1.0,true,b,true,"pfizer",30,false)#15 extra doses, keeping priority groups
    run_param_fix([20],1.0,1.0,true,b,true,"pfizer",30,true)#15 extra doses given to general population
    run_param_fix([20],1.0,1.0,true,b,true,"pfizer",45,false)#15 extra doses, keeping priority groups
    run_param_fix([20],1.0,1.0,true,b,true,"pfizer",45,true)#15 extra doses given to general population
 
    run_param_fix([20],1.0,1.0,true,b,true,"pfizer",60,false)#15 extra doses, keeping priority groups
    run_param_fix([20],1.0,1.0,true,b,true,"pfizer",60,true)#15 extra doses given to general population

    ###Moderna
    #no reduction
    
    #50% reduction
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",0,false)#no extra dose at day 80
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",15,false)#15 extra doses, keeping priority groups
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",15,true)#15 extra doses given to general population
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",30,false)#15 extra doses, keeping priority groups
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",30,true)#15 extra doses given to general population
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",45,false)#15 extra doses, keeping priority groups
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",45,true)#15 extra doses given to general population 
    
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",60,false)#15 extra doses, keeping priority groups
    run_param_fix([20],1.0,1.0,true,b,true,"moderna",60,true)#15 extra doses given to general population

end 