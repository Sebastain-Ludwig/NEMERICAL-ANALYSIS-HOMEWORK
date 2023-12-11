function root=nnl_dfr(left_point,right_point,function_dfr,tor_err)
%%%DIVIDE FIND ROOT OF A FUNCTION
%%%FURE FUNCTION

%%  FUNCTION_DFR
    runsimulation=function_dfr;

%%  SECTION:SPAN INITIATE
    lp=left_point;%SPAN:[lp,rp]
    rp=right_point;%SPAN:[lp,rp]
    lpv=runsimulation(left_point);%FUN VAR:lpv=f(lp)
    rpv=runsimulation(right_point);%FUN VAR:rpv:f(rp)
    
%%  SECTION:COMPLEXITY ESTIMATE 
    sim_n=ceil(log2((rp-lp)/tor_err));%FIND STEP NUMBER CALCULATE
    disp(["THE WORST SITUATION ITER STEP:",sim_n,"SIMULATION NUMBER:",sim_n+1]);
 
%%  SECTION:DIVIDE FIND ROOT
    mp=(lp+rp)/2.0;%MID POINT SPAN:[lp,mp=(lp+rp)/2.0,rp]
    mdv=runsimulation((lp+rp)/2.0);%TIME CONSUMING:<1>
    mdv_vec=[];
    
    notrepeat=@(x,vec)not(ismember(x,vec(1:end-1)));
    
    while(abs(mdv)>tor_err&&notrepeat(mdv,mdv_vec))%BREAK WHIN <TOR_ERROR> IS SATISFIED
        mp=(lp+rp)/2.0;
        mdv=runsimulation(mp);%TIME CONSUMING:<1>
        mdv_vec=[mdv_vec,mdv];
        
        if(lpv*mdv<0)%f(lp)*f(mp)? 0
            rp=mp;
            rpv=mdv;
        else
            lp=mp;
            lpv=mdv;
        end
    end
    root=mp;
end