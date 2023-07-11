function get_beta(t)
    ind = 1;
    
    # first day of simulation
    # day0 = "1-jan-2020";

    tlock = zeros(13)
    tunlock = zeros(13)
    scale = zeros(13)
    
    # state of emergency and such [2]
    tlock[ind]   = 82;  # daysact(day0,'23-mar-2020');
    tunlock[ind] = 198; # daysact(day0,'17-jul-2020');
    scale[ind]   = 0.4;
    ind = ind+1;
    
    # stage 3 [3]
    tlock[ind]   = 198; # daysact(day0,'17-jul-2020');
    tunlock[ind] = 282; # daysact(day0,'9-oct-2020');
    scale[ind]   = 0.55;
    ind = ind+1;
    
    # oops, back to stage 2 [2]
    tlock[ind]   = 282; # daysact(day0,'9-oct-2020');
    tunlock[ind] = 327; # daysact(day0,'23-nov-2020');
    scale[ind]   = 0.6;
    ind = ind+1;
    
    # strict lockdown, stay at home order [1]
    tlock[ind]   = 327; # daysact(day0,'23-nov-2020');
    # tunlock[ind] = daysact(day0,'5-mar-2021');
    tunlock[ind] = 370; # daysact(day0,'5-jan-2021');
    scale[ind]   = 0.55;
    ind = ind+1;
    
    # stricter lockdown, stay at home order [1]
    tlock[ind]   = 370; # daysact(day0,'5-jan-2021');
    tunlock[ind] = 429; # daysact(day0,'5-mar-2021');
    scale[ind]   = 0.3; # 0.25;
    ind = ind+1;
    
    # stay at home order lifted [2]
    tlock[ind]   = 429; # daysact(day0,'5-mar-2021');
    tunlock[ind] = 468; # daysact(day0,'13-apr-2021');
    scale[ind]   = 0.4; # 0.4;
    ind = ind+1;
    
    # lockdown again [1]
    tlock[ind]   = 468; # daysact(day0,'13-apr-2021');
    tunlock[ind] = 527; # daysact(day0,'11-jun-2021');
    scale[ind]   = 0.15; # 0.1;
    ind = ind+1;
    
    # stage 2 reopening [3]
    tlock[ind]   = 527; #daysact(day0,'11-jun-2021');
    tunlock[ind] = 548; #daysact(day0,'2-jul-2021');
    scale[ind]   = 0.2; #0.15;
    ind = ind+1;
    
    # stage 3 reopening [4]
    tlock[ind]   = 548; # daysact(day0,'2-jul-2021');
    tunlock[ind] = 562; # daysact(day0,'16-jul-2021');
    scale[ind]   = 0.3; # 0.25;
    ind = ind+1;
    
    # physical distancing [5]
    tlock[ind]   = 562; # daysact(day0,'16-jul-2021');
    tunlock[ind] = 609; # daysact(day0,'1-sep-2021');
    scale[ind]   = 0.35;
    ind = ind+1;
    
    # physical distancing [5]
    tlock[ind]   = 609; # daysact(day0,'1-sep-2021');
    tunlock[ind] = 670; # daysact(day0,'1-nov-2021');
    scale[ind]   = 0.4;
    ind = ind+1;
    
    # physical distancing loosening [6]
    tlock[ind]   = 670; # daysact(day0,'1-nov-2021');
    tunlock[ind] = 731; # daysact(day0,'1-jan-2022');
    scale[ind]   = 0.5;
    ind = ind+1;
    
    # physical distancing tightened [6]
    tlock[ind]   = 731; # daysact(day0,'1-jan-2022');
    tunlock[ind] = 1096; # daysact(day0,'1-jan-2023');
    scale[ind]   = 0.4;
    ind = ind+1;
    
    ilock  = findfirst(tlock.<=t);
    iunlock= findfirst(t.<tunlock);

    @show t,tlock,ilock,iunlock,tunlock

    beta_scale = if ilock !== nothing && iunlock !== nothing
        ii = intersect(ilock,iunlock);
        if ii === nothing
            1
        else
            scale[ii]
        end
    else
        1
    end
end

function svair_simple(dy,y,p,t)
    
    beta,beta_v1,beta_v2,beta_R,ai_beta_ratio,gamma,nu_v1,nu_v2,nu_R,ai,ai_V,ai_R,mu,mu_I,mu_IV = p

    # retrieve current populations
    n_var = 1;  # number of viruses simulated
    n_vax = 1;
    ind  = 1;
    S    = y[ind]; ind=ind+1;  # original susceptible
    SVR  = y[ind]; ind=ind+1;  # lost immunity after vaccination or recovery
    V1   = y[ind:ind+n_vax-1]; ind=ind+n_vax;  # one-dose vaccination
    V2   = y[ind:ind+n_vax-1]; ind=ind+n_vax;  # fully vaccinated
    I    = y[ind:ind+n_var-1]; ind=ind+n_var;  # infected
    IV   = y[ind:ind+n_var-1]; ind=ind+n_var;  # infected even with vaccination
    IR   = y[ind:ind+n_var-1]; ind=ind+n_var;  # infected again after recovery from a different variant
    A    = y[ind:ind+n_var-1]; ind=ind+n_var;  # asymptomatic infections
    AR   = y[ind:ind+n_var-1]; ind=ind+n_var;  # asymptomatic infections after recovery from a different variant
    R    = y[ind:ind+n_var-1]; ind=ind+n_var;  # recovered
    R2   = y[ind]; ind=ind+1;  # recovered after getting both variants

    # get time-dependent parameters
    # [vr1, vr2] = get_vaccine_rate (t);
    # no vaccine for the first year
    vr1 = 0;
    vr2 = 0;

    beta_scale = get_beta(t);
    beta       = beta*beta_scale;
    beta_v1    = beta_v1*beta_scale;
    beta_v2    = beta_v2*beta_scale;
    beta_R     = beta_R*beta_scale;

    @show I,IV,IR,ai_beta_ratio
    @show I+IV+IR+ai_beta_ratio
    @show A+AR
    # total infectious population
    I_total = I+IV+IR+ai_beta_ratio.*(A+AR);

    ## need the following to compute infection of recovered from another variant
    # mm = ones(n_var+1)-diag(ones(1,n_var+1));
    # mv = mm.*repmat(R,1,n_var+1);
    Rv = 0; #sum(mv)'; # in simple case, don't consider infection from a different variant

    @show vr1.*(S+sum(A)) 
    @show vr2.*V1 
    @show nu_v1.*V1 
    @show sum(beta_v1.*(V1*(I_total)'),dims=2) 
    @show mu*V1

    # compute time derivatives
    dSdt   = - sum(beta.*S.*(I_total)) - sum(vr1.*S) + mu*(1-S);
    dSVRdt = + nu_v1.*sum(V1) + nu_v2.*sum(V2) + sum(nu_R.*R) + nu_R.*R2 - sum(beta.*SVR.*(I_total)) - mu*SVR;
    dV1dt  = + vr1.*(S+sum(A)) .- vr2.*V1 .- nu_v1.*V1 .- sum(beta_v1.*(V1*(I_total)'),dims=2) .- mu*V1;
    dV2dt  = + vr2.*V1 - nu_v2.*V2 - sum(beta_v2.*(V2.*(I_total)'),dims=2) - mu*V2;
    dIdt   = (1-ai).*(+ beta.*S.*(I_total) + beta.*SVR.*(I_total)) - gamma.*I - mu_I.*I;
    dIVdt  = (1-ai_V).*(+ sum(beta_v1.*(V1*(I_total)'))' + sum(beta_v2.*(V2.*(I_total)'))') - gamma.*IV - mu_IV.*IV;
    dIRdt  = (1-ai_R).*(+ beta_R.*Rv.*(I_total)) - gamma.*IR - mu*IR;
    dAdt   = ai.*(+ beta.*S.*(I_total) + beta.*SVR.*(I_total)) + ai_V.*(sum(beta_v1.*(V1*(I_total)'))' + sum(beta_v2.*(V2.*(I_total)'))') - sum(vr1).*A - gamma.*A - mu.*A;
    dARdt  = ai_R.*(+ beta_R.*Rv.*(I_total)) - gamma.*AR - mu*AR;
    dRdt   = + gamma.*(I_total) - nu_R.*R - beta_R.*Rv.*(I_total) - mu*R;
    dR2dt  = + sum(gamma.*(IR+AR)) - nu_R.*R2 - mu*R2;

    dy .= [dSdt;dSVRdt;dV1dt;dV2dt;dIdt;dIVdt;dIRdt;dAdt;dARdt;dRdt;dR2dt];
end

# define integration interval
t0 = 0;
tfinal = 2.0*365;

# Ontario population, 14.57M 2019
N = 14570000;

# define parameters
# these are infection rates without lockdown
beta = [3.3e-9;5.5e-9;7.6e-9]*N; # infection rate, for susceptibles
# beta_v1 = [0.3;0.5].*beta; % infection rate, first dose
# beta_v2 = 0.05*beta; % infection rate, both doses
#% https://www.gov.uk/government/news/vaccines-highly-effective-against-b-1-617-2-variant-after-2-doses
beta_v1 = [0.2 0.5 0.67; 0.2 0.5 0.67].*[beta';beta']; # infection rate, first dose; row, for a given vaccine type; column, for a given variant
beta_v2 = [0.05 0.07 0.12; 0.05 0.34 0.4].*[beta';beta']; # infection rate, both doses
beta_R  = 0.05*beta; # infection rate after recovery
ai_beta_ratio = [3] # [3; 3; 3];  # asymptomatic vs. symptomatic infectivity ratios
# vaccination rates now given as function of time in get_vaccine_rates.m
# vr1  = 1e-3; % vaccination rates (per day)
# vr2  = 1/21; % 21 days delay
gamma = 1/28; # recovery rate
nu_v1 = 2*0.25/182; # loss of immunity, first dose (6 months)
nu_v2 = 2*0.125/365; # loss of immunity, both doses (1 year)
nu_R  = 2*0.125/365; # loss of immunity, recovered (1 year)
ai    = [0.5; 0.5; 0.5];  # fraction of asymtomatic primary infections
ai_V  = [0.85; 0.85; 0.85];  # fraction of asymtomatic infections after vacciation
ai_R  = [0.85; 0.85; 0.85];  # fraction of asymtomatic infections after recovery from another variant
mu    = 109019/N/365;  # natural death rate (109019 in 14.5 M in 2018-2019)
mu_I  = 1.75*[9255/555927*gamma; 1.6*9255/555927*gamma; 1.8*9255/555927*gamma];  # COVID mortaolity rate, (9255 deaths for 555927 total cases)
mu_IV = 0.15*mu_I;  # vaccine reduces mortality rate

p = (beta,beta_v1,beta_v2,beta_R,ai_beta_ratio,gamma,nu_v1,nu_v2,nu_R,ai,ai_V,ai_R,mu,mu_I,mu_IV)

# parameters for new killer variant, will replace wild-type after fall 2021
new_beta    = beta[3]; #2.2*beta(1);
new_beta_v1 = [0.5; 0.5]*new_beta;
new_beta_v2 = [0.2; 0.2]*new_beta;
new_beta_R  = 0.05*new_beta;
new_ai      = 0.8;
t_new_voc   = daysact('1-jan-2020','1-sep-2022');

# define initial population fractions
I0   = [1e-6;0;0];  # infected
A0   = [0;0;0];
S0   = 1-sum(I0+A0);  # original susceptible
SVR0 = 0;  # lost immunity after vaccination or recovery
V10  = [0;0];  # one-dose vaccination
V20  = [0;0];  # fully vaccinated
IV0  = [0;0;0];  # infected even with vaccination
IR0  = [0;0;0];  # infected again after recovery from a different variant
AR0  = [0;0;0];  # asymptomatic infection after recovery from a different variant
R0   = [0;0;0];  # recovered
R20  = 0;  # recovered after getting both variants
y0 = [S0;SVR0;V10;V20;I0;IV0;IR0;A0;AR0;R0;R20];

using OrdinaryDiffEq
prob = ODEProblem(svair_simple, y0, (t0,tfinal), p)
sol = solve(prob, Tsit5())
