## Core equations for transitions between compartments:
# Notation: XY where X=PWID status[A(active), E(ex will relapse), Q(ex will not relapse)], Y=HCV status[N(negative), A(acute), C(chronic)]
AN_ = AN + n_debut_AN - n_AN_EN - n_AN_AA + n_EN_AN - n_AN_D - n_AN_QN #
AA_ = AA + n_AN_AA - n_AA_AR - n_AA_AC - n_AA_EA + n_EA_AA - n_AA_QA + n_AR_AA - n_AA_D             
AC_ = AC + n_AA_AC - n_AC_EC + n_EC_AC - n_AC_AR - n_AC_D - n_AC_QC + n_IC_AC           
AR_ = AR + n_AA_AR + n_AC_AR - n_AR_ER + n_ER_AR - n_AR_QR - n_AR_AA - n_AR_D                                           
EN_ = EN + n_AN_EN - n_EN_AN - n_EN_D                         
EA_ = EA + n_AA_EA - n_EA_EC - n_EA_AA - n_EA_ER - n_EA_D                         
EC_ = EC + n_AC_EC + n_EA_EC - n_EC_ER - n_EC_AC - n_EC_D                       
ER_ = ER + n_EC_ER + n_EA_ER + n_AR_ER - n_ER_AR - n_ER_D                                           
QN_ = QN + n_AN_QN - n_QN_D                         
QA_ = QA + n_AA_QA - n_QA_QC - n_QA_QR - n_QA_D                        
QC_ = QC + n_AC_QC + n_QA_QC - n_QC_QR - n_QC_D                      
QR_ = QR + n_QC_QR + n_QA_QR + n_AR_QR - n_QR_D                                            
D_ = D + n_AN_D + n_AA_D + n_AC_D + n_AR_D + n_EN_D + n_EA_D + n_EC_D + n_ER_D + n_QN_D + n_QA_D + n_QC_D + n_QR_D
# Immigrant compartments:
IC_ = IC + n_immigrated_hcvpos - n_IC_IR - n_IC_D - n_IC_AC
IR_ = IR + n_IC_IR - n_IR_D
ID_ = ID + n_IC_D + n_IR_D
# General population compartments
GC_ = GC + n_genpop_hcvpos - n_GC_GR - n_GC_D
GR_ = GR + n_GC_GR - n_GR_D
GD_ = GD + n_GC_D
# Use threshold
upper_threshold = 2**25 #  -- upper threshold to avoid overflow 
AN__ = if(AN_ < upper_threshold) AN_ else upper_threshold
AA__ = if(AA_ < upper_threshold) AA_ else upper_threshold
AC__ = if(AC_ < upper_threshold) AC_ else upper_threshold
AR__ = if(AR_ < upper_threshold) AR_ else upper_threshold
EN__ = if(EN_ < upper_threshold) EN_ else upper_threshold
EA__ = if(EA_ < upper_threshold) EA_ else upper_threshold
EC__ = if(EC_ < upper_threshold) EC_ else upper_threshold
ER__ = if(ER_ < upper_threshold) ER_ else upper_threshold
QN__ = if(QN_ < upper_threshold) QN_ else upper_threshold
QA__ = if(QA_ < upper_threshold) QA_ else upper_threshold
QC__ = if(QC_ < upper_threshold) QC_ else upper_threshold
QR__ = if(QR_ < upper_threshold) QR_ else upper_threshold
D__ = if(D_ < upper_threshold) D_ else upper_threshold
IC__ = if(IC_ < upper_threshold) IC_ else upper_threshold
IR__ = if(IR_ < upper_threshold) IR_ else upper_threshold
ID__ = if(ID_ < upper_threshold) ID_ else upper_threshold
GC__ = if(GC_ < upper_threshold) GC_ else upper_threshold
GR__ = if(GR_ < upper_threshold) GR_ else upper_threshold
GD__ = if(GD_ < upper_threshold) GD_ else upper_threshold

# Non-negative threshold:
update(AN) = if(AN__ > 0) AN__ else 0
update(AA) = if(AA__ > 0) AA__ else 0
update(AC) = if(AC__ > 0) AC__ else 0
update(AR) = if(AR__ > 0) AR__ else 0
update(EN) = if(EN__ > 0) EN__ else 0
update(EA) = if(EA__ > 0) EA__ else 0
update(EC) = if(EC__ > 0) EC__ else 0
update(ER) = if(ER__ > 0) ER__ else 0
update(QN) = if(QN__ > 0) QN__ else 0
update(QA) = if(QA__ > 0) QA__ else 0
update(QC) = if(QC__ > 0) QC__ else 0
update(QR) = if(QR__ > 0) QR__ else 0
update(D) = if(D__ > 0) D__ else 0
update(IC) = if(IC__ > 0) IC__ else 0
update(IR) = if(IR__ > 0) IR__ else 0
update(ID) = if(ID__ > 0) ID__ else 0
update(GC) = if(GC__ > 0) GC__ else 0
update(GR) = if(GR__ > 0) GR__ else 0
update(GD) = if(GD__ > 0) GD__ else 0

# Random walk on rates to be inferred by particle filter
rd_multiplier_ = rnorm(1, sd_rw_debut)
rd_multiplier__ = if (rd_multiplier_ > 1e-4) rd_multiplier_ else 1e-4
rd_multiplier = if (((step+1)*dt+1) >= time_forecast) 1 else rd_multiplier__ # No random walk after time_forecast
update(rate_debut) = rate_debut * rd_multiplier
# Rate of active PWD quitting to will-not-relapse state (all AY -> QY)
update(rate_quitting_wnr) = rate_quitting_wnr
# Rate of active PWD quitting to will-relapse state (all AY -> EY)
update(rate_quitting_wr) = rate_quitting_wr
# Rate of ex PWID relapsing (all EY -> AY)
update(rate_relapsing) = rate_relapsing
# Beta transmission parameter for force of infection, causing AN -> AA
update(beta) = beta 
# Rate of recovery
update(rate_recovery_ex) = rate_recovery_ex # exp(log_update_rrecex)
# Rate of treatment
rtreat_multiplier_ = rnorm(1, sd_rtreat)
rtreat_multiplier__ = if (rtreat_multiplier_ > 1e-4) rtreat_multiplier_ else 1e-4
rtreat_multiplier = if (((step+1)*dt+1) >= time_forecast) 1 else rtreat_multiplier__ # No random walk after time_forecast
# We want to impose a bounded random walk on the rate_treatment parameter, but we need to do it in the unbounded space, so we need to transform it to the unbounded space, do the random walk, and then transform it back to the bounded space
rtreat_unbounded = log(rate_treatment / (rate_treatment_max - rate_treatment)) # Logit transformation to unbounded space
rtreat_unbounded_update = rtreat_unbounded * rtreat_multiplier # Multiplicative random walk in unbounded space
rtreat_bounded_update = rate_treatment_max / (1 + exp(-rtreat_unbounded_update)) # Transform back to bounded space
# The complicated piecewise thing for the rate of treament, which suddenly jumps to the rate_treatment_rw_ini value at time_start_rw_treatment, then random walks from there
update(rate_treatment) = if (switched_to_rw_treatment < 1) rate_treatment_input_currstep else if (switched_to_rw_treatment == 1) rate_treatment_rw_ini else rtreat_bounded_update
update(switched_to_rw_treatment) = if (((step+1) * dt + 1) >= time_start_rw_treatment) switched_to_rw_treatment+1 else 0
update(rate_treatment_active) = rate_treatment_active_input_currstep


## Individual probabilities of transition:
p_quitting_wnr = 1 - exp(-rate_quitting_wnr * dt) # will not relapse
p_quitting_wr = 1 - exp(-rate_quitting_wr * dt) # will relapse
p_relapsing = 1 - exp(-rate_relapsing * dt)
p_FOI_ = 1 - exp(- beta * ((if ((AA+AC) > 0) (AA+AC) else 0) / (if (N_active > 1) N_active else 1)) * dt * gini_drug_deaths) # Force of infection, following Meijerink et al 
p_FOI = p_FOI_ * risk_multiplier_lar_nsp # Combined effect of LAR and NSP coverage, reducing FOI. (1 - nsp_coverage * max_rel_reduction_nsp) # Effect of NSP modifying FOI, following Meijerink et al
p_recover_ex = 1 - exp(-rate_recovery_ex * dt) # Prob recovery ex-PWID
p_recover_active = p_recover_ex * rel_prob_recovery_active # Prob recovery active PWID
p_treatment_ex = 1 - exp(-rate_treatment * dt) 
p_treatment_active = 1 - exp(-rate_treatment_active * dt)
p_acute_chronic = 1 - exp(-rate_acute_chronic * dt)
p_death_ex = 1 - exp(-rate_death_ex * dt)
p_death_active = 1 - exp(-rate_death_active * dt)
# Immigrants and gen pop
p_treatment_imm = 1 - exp(-rate_treatment * dt) 
p_death_imm = 1 - exp(-rate_death_imm * dt)
p_treatment_genpop = 1 - exp(-rate_treatment * dt) 
p_death_genpop = 1 - exp(-rate_death_genpop * dt)
p_IC_AC = 1 - exp(-rate_IC_AC * dt)

## Draws from binomial distributions for numbers changing between compartments:
n_debut_AN_ = rpois(rate_debut*dt)                                            # Debut of new individuals as active PWID, AN. Poisson because the rate does not depend on the number of active PWID.
n_AN_EN_ = rbinom(AN + n_debut_AN_, p_quitting_wr)                            # Active HCV negative PWID quitting
n_EN_AN_ = rbinom(EN + n_AN_EN_, p_relapsing)                              # Ex-PWID, HCV negative, relapsing
n_AN_AA_ = rbinom(AN + n_debut_AN_ - n_AN_EN_ + n_EN_AN_, p_FOI)                       # Active PWID, HCV naive, getting infected by HCV (acute)
n_AR_AA_ = rbinom(AR, rel_prob_reinfection * p_FOI)                                                  # Active PWID, recovered, getting reinfected by HCV (acute)
n_AA_EA_ = rbinom(AA + n_AN_AA_ + n_AR_AA_, p_quitting_wr)                               # Active PWID, acutely infected, quitting injection
n_AA_AR_ = rbinom(AA + n_AN_AA_ + n_AR_AA_ - n_AA_EA_, p_recover_active)               # Active PWID, acutely infected, spontaneously clearing infection
n_AA_AC_ = rbinom(AA + n_AN_AA_ + n_AR_AA_ - n_AA_EA_ - n_AA_AR_, p_acute_chronic)      # Active PWID, acute infection becoming chronic
n_EA_AA_ = rbinom(EA + n_AA_EA_, p_relapsing)                              # Ex-PWID, acutely infected, relapsing to active PWID  
n_EA_EC_ = rbinom(EA + n_AA_EA_ - n_EA_AA_, p_acute_chronic)                # Ex-PWID, acute infection becoming chronic
n_AC_EC_ = rbinom(AC + n_AA_AC_, p_quitting_wr)                               # Active PWID, chronically infected, quitting PWID
n_AC_AR_ = rbinom(AC + n_AA_AC_ - n_AC_EC_, p_treatment_active*p_treatment_success)  # Active PWID, chronically infected, getting treated to clear infection
n_EC_ER_ = rbinom(EC + n_EA_EC_ + n_AC_EC_, p_treatment_ex*p_treatment_success)  # Ex-PWID, chronically infected, getting treated to clear infection
n_EC_AC_ = rbinom(EC + n_EA_EC_ + n_AC_EC_ - n_EC_ER_, p_relapsing)          # Ex-PWID, chronically infected, relapsing to become active PWID
n_EA_ER_ = rbinom(EA + n_AA_EA_ - n_EA_AA_ - n_EA_EC_, p_recover_ex)         # Active PWID, acutely infected, spontaneously clearing infection
n_AN_QN_ = rbinom(AN + n_debut_AN_ - n_AN_EN_ + n_EN_AN_ - n_AN_AA_, p_quitting_wnr) # Active PWID, HCV negative, quitting without relapsing
n_QC_QR_ = rbinom(QC, p_treatment_ex*p_treatment_success)  # Quitters, chronically infected, getting treated to clear infection
n_QA_QR_ = rbinom(QA, p_recover_ex)  # Quitters, acutely infected, getting treated to clear infection
n_AA_QA_ = rbinom(AA + n_AN_AA_ + n_AR_AA_ - n_AA_EA_ - n_AA_AR_ - n_AA_AC_ + n_EA_AA_, p_quitting_wnr) # Active PWID, acutely infected, quitting injection permanently
n_QA_QC_ = rbinom(QA - n_QA_QR_ + n_AA_QA_, p_acute_chronic) # Quitters, acutely infected, acute infection becoming chronic
n_AC_QC_ = rbinom(AC + n_AA_AC_ - n_AC_EC_ - n_AC_AR_ + n_EC_AC_, p_quitting_wnr) # Active PWID, chronically infected, quitting PWID and relapsing
n_AR_ER_ = rbinom(AR - n_AR_AA_ + n_AA_AR_ + n_AC_AR_, p_quitting_wr) # Active PWID, recovered, quitting to wr
n_ER_AR_ = rbinom(ER + n_EC_ER_ + n_EA_ER_ + n_AR_ER_, p_relapsing) # Ex-PWID, recovered, relapsing to active PWID
n_AR_QR_ = rbinom(AR -n_AR_AA_ + n_AA_AR_ + n_AC_AR_ - n_AR_ER_, p_quitting_wnr) # Active PWID, recovered, quitting to wnr
n_AN_D_ = rbinom(AN + n_debut_AN_ - n_AN_EN_ + n_EN_AN_ - n_AN_AA_ - n_AN_QN_, p_death_active) # Active PWID, HCV negative, dying
n_EN_D_ = rbinom(EN + n_AN_EN_ - n_EN_AN_, p_death_ex) # Ex-PWID, HCV negative, dying
n_AA_D_ = rbinom(AA + n_AN_AA_ + n_AR_AA_ - n_AA_EA_ - n_AA_AR_ - n_AA_AC_ + n_EA_AA_ - n_AA_QA_, p_death_active) # Active PWID, acutely infected, dying
n_AC_D_ = rbinom(AC + n_AA_AC_ - n_AC_EC_ - n_AC_AR_ + n_EC_AC_ - n_AC_QC_, p_death_active) # Active PWID, chronically infected, dying
n_AR_D_ = rbinom(AR - n_AR_AA_ + n_AA_AR_ + n_AC_AR_ - n_AR_ER_ + n_ER_AR_ - n_AR_QR_, p_death_active) # Active PWID, recovered, dying
n_EA_D_ = rbinom(EA + n_AA_EA_ - n_EA_AA_ - n_EA_EC_ - n_EA_ER_, p_death_ex) # Ex-PWID, acutely infected, dying
n_EC_D_ = rbinom(EC + n_EA_EC_ + n_AC_EC_ - n_EC_ER_ - n_EC_AC_, p_death_ex) # Ex-PWID, chronically infected, dying
n_ER_D_ = rbinom(ER + n_EC_ER_ + n_EA_ER_ + n_AR_ER_ - n_ER_AR_, p_death_ex) # Ex-PWID, recovered, dying
n_QN_D_ = rbinom(QN + n_AN_QN_, p_death_ex) # Quitters, HCV negative, dying
n_QA_D_ = rbinom(QA - n_QA_QR_ + n_AA_QA_ - n_QA_QC_, p_death_ex) # Quitters, acutely infected, dying
n_QC_D_ = rbinom(QC - n_QC_QR_ + n_QA_QC_ + n_AC_QC_, p_death_ex) # Quitters, chronically infected, dying
n_QR_D_ = rbinom(QR + n_QC_QR_ + n_QA_QR_ + n_AR_QR_, p_death_ex) # Quitters, recovered, dying
# Immigrants:
n_immigrated_hcvpos = immigrated_hcvpos_currstep
n_IC_IR_ = rbinom(IC + n_immigrated_hcvpos, p_treatment_imm*p_treatment_success)
n_IC_AC_ = rbinom(IC + n_immigrated_hcvpos - n_IC_IR_, p_IC_AC) # Immigrants debuting as PWID, chronically infected
n_IC_D_ = rbinom(IC + n_immigrated_hcvpos - n_IC_IR_ - n_IC_AC_, p_death_imm)
n_IR_D_ = rbinom(IR + n_IC_IR_, p_death_imm)
# General population:
n_genpop_hcvpos = genpop_hcvpos_currstep
n_GC_GR_ = rbinom(GC + n_genpop_hcvpos, p_treatment_genpop*p_treatment_success)
n_GC_D_ = rbinom(GC + n_genpop_hcvpos - n_GC_GR_, p_death_genpop)
n_GR_D_ = rbinom(GR + n_GC_GR_, p_death_genpop)
# Checks to avoid overflow: (rather do it here than directly on compartments above, to conserve population size)
upper_threshold_trans = 1e3 
n_debut_AN = if (n_debut_AN_ < upper_threshold_trans) n_debut_AN_ else upper_threshold_trans
n_AN_EN = if (n_AN_EN_ < upper_threshold_trans) n_AN_EN_ else upper_threshold_trans
n_EN_AN = if (n_EN_AN_ < upper_threshold_trans) n_EN_AN_ else upper_threshold_trans
n_AN_AA = if (n_AN_AA_ < upper_threshold_trans) n_AN_AA_ else upper_threshold_trans
n_AR_AA = if (n_AR_AA_ < upper_threshold_trans) n_AR_AA_ else upper_threshold_trans
n_AA_EA = if (n_AA_EA_ < upper_threshold_trans) n_AA_EA_ else upper_threshold_trans
n_AA_AR = if (n_AA_AR_ < upper_threshold_trans) n_AA_AR_ else upper_threshold_trans
n_AA_AC = if (n_AA_AC_ < upper_threshold_trans) n_AA_AC_ else upper_threshold_trans
n_EA_AA = if (n_EA_AA_ < upper_threshold_trans) n_EA_AA_ else upper_threshold_trans
n_EA_EC = if (n_EA_EC_ < upper_threshold_trans) n_EA_EC_ else upper_threshold_trans
n_AC_EC = if (n_AC_EC_ < upper_threshold_trans) n_AC_EC_ else upper_threshold_trans
n_AC_AR = if (n_AC_AR_ < upper_threshold_trans) n_AC_AR_ else upper_threshold_trans
n_EC_ER = if (n_EC_ER_ < upper_threshold_trans) n_EC_ER_ else upper_threshold_trans
n_EC_AC = if (n_EC_AC_ < upper_threshold_trans) n_EC_AC_ else upper_threshold_trans
n_EA_ER = if (n_EA_ER_ < upper_threshold_trans) n_EA_ER_ else upper_threshold_trans
n_AN_QN = if (n_AN_QN_ < upper_threshold_trans) n_AN_QN_ else upper_threshold_trans
n_QC_QR = if (n_QC_QR_ < upper_threshold_trans) n_QC_QR_ else upper_threshold_trans
n_QA_QR = if (n_QA_QR_ < upper_threshold_trans) n_QA_QR_ else upper_threshold_trans
n_AA_QA = if (n_AA_QA_ < upper_threshold_trans) n_AA_QA_ else upper_threshold_trans
n_QA_QC = if (n_QA_QC_ < upper_threshold_trans) n_QA_QC_ else upper_threshold_trans
n_AC_QC = if (n_AC_QC_ < upper_threshold_trans) n_AC_QC_ else upper_threshold_trans
n_AR_ER = if (n_AR_ER_ < upper_threshold_trans) n_AR_ER_ else upper_threshold_trans
n_ER_AR = if (n_ER_AR_ < upper_threshold_trans) n_ER_AR_ else upper_threshold_trans
n_AR_QR = if (n_AR_QR_ < upper_threshold_trans) n_AR_QR_ else upper_threshold_trans
n_AN_D = if (n_AN_D_ < upper_threshold_trans) n_AN_D_ else upper_threshold_trans
n_AA_D = if (n_AA_D_ < upper_threshold_trans) n_AA_D_ else upper_threshold_trans
n_AC_D = if (n_AC_D_ < upper_threshold_trans) n_AC_D_ else upper_threshold_trans
n_AR_D = if (n_AR_D_ < upper_threshold_trans) n_AR_D_ else upper_threshold_trans
n_EN_D = if (n_EN_D_ < upper_threshold_trans) n_EN_D_ else upper_threshold_trans
n_EA_D = if (n_EA_D_ < upper_threshold_trans) n_EA_D_ else upper_threshold_trans
n_EC_D = if (n_EC_D_ < upper_threshold_trans) n_EC_D_ else upper_threshold_trans
n_ER_D = if (n_ER_D_ < upper_threshold_trans) n_ER_D_ else upper_threshold_trans
n_QN_D = if (n_QN_D_ < upper_threshold_trans) n_QN_D_ else upper_threshold_trans
n_QA_D = if (n_QA_D_ < upper_threshold_trans) n_QA_D_ else upper_threshold_trans
n_QC_D = if (n_QC_D_ < upper_threshold_trans) n_QC_D_ else upper_threshold_trans
n_QR_D = if (n_QR_D_ < upper_threshold_trans) n_QR_D_ else upper_threshold_trans
# Immigrants:
n_IC_IR = if (n_IC_IR_ < upper_threshold_trans) n_IC_IR_ else upper_threshold_trans
n_IC_AC = if (n_IC_AC_ < upper_threshold_trans) n_IC_AC_ else upper_threshold_trans
n_IC_D = if (n_IC_D_ < upper_threshold_trans) n_IC_D_ else upper_threshold_trans
n_IR_D = if (n_IR_D_ < upper_threshold_trans) n_IR_D_ else upper_threshold_trans
# General population:
n_GC_GR = if (n_GC_GR_ < upper_threshold_trans) n_GC_GR_ else upper_threshold_trans
n_GC_D = if (n_GC_D_ < upper_threshold_trans) n_GC_D_ else upper_threshold_trans
n_GR_D = if (n_GR_D_ < upper_threshold_trans) n_GR_D_ else upper_threshold_trans

## Update compartment counts:
update(cum_infections) = cum_infections + n_AN_AA
update(inc_infections) = if (step %% as.integer(1/dt) == 0) n_AN_AA else inc_infections + n_AN_AA # Only update incidence on integer timesteps, cumulatively sum between timesteps
update(cum_treatments) = cum_treatments + n_AC_AR + n_EC_ER
update(inc_treatments_pwid) = if (step %% as.integer(1/dt) == 0) as.integer((n_AC_AR + n_EC_ER + n_QC_QR)/p_treatment_success) else inc_treatments_pwid + as.integer((n_AC_AR + n_EC_ER + n_QC_QR)/p_treatment_success)
update(cum_debut) = cum_debut + n_debut_AN
update(inc_debut) = if (step %% as.integer(1/dt) == 0) n_debut_AN else inc_debut + n_debut_AN
update(inc_death_active) = if (step %% as.integer(1/dt) == 0) n_AN_D + n_AA_D + n_AC_D else inc_death_active + n_AN_D + n_AA_D + n_AC_D
update(inc_death_ex_wr) = if (step %% as.integer(1/dt) == 0) n_EN_D + n_EA_D + n_EC_D else inc_death_ex_wr + n_EN_D + n_EA_D + n_EC_D
update(inc_death_ex_wnr) = if (step %% as.integer(1/dt) == 0) n_QN_D + n_QA_D + n_QC_D else inc_death_ex_wnr + n_QN_D + n_QA_D + n_QC_D
update(inc_death_immigrants) = if (step %% as.integer(1/dt) == 0) n_IC_D else inc_death_immigrants + n_IC_D
update(inc_death_genpop) = if (step %% as.integer(1/dt) == 0) n_GC_D else inc_death_genpop + n_GC_D
update(inc_treatments_immigrants) = if (step %% as.integer(1/dt) == 0) as.integer((n_IC_IR)/p_treatment_success) else inc_treatments_immigrants + as.integer((n_IC_IR)/p_treatment_success)
update(inc_treatments_genpop) = if (step %% as.integer(1/dt) == 0) as.integer((n_GC_GR)/p_treatment_success) else inc_treatments_genpop + as.integer((n_GC_GR)/p_treatment_success)
update(inc_immigrated_hcvpos) = if (step %% as.integer(1/dt) == 0) n_immigrated_hcvpos else inc_immigrated_hcvpos + n_immigrated_hcvpos # TODO understand if this is an off-by-one-timestep error -- should it be (step+1) %% as.integer(1/dt) == 0? Except that produces insane results with much lower incidence counts.
update(cum_immigrated_hcvpos) = cum_immigrated_hcvpos + n_immigrated_hcvpos
update(cum_IC_AC) = cum_IC_AC + n_IC_AC
update(cum_AA_AC) = cum_AA_AC + n_AA_AC

# Update the interpolation vector using workaround from https://mrc-ide.github.io/odin.dust/articles/porting.html
update(p_treatment_success) = if (as.integer(step) >= length(p_treatment_succcess_step)) p_treatment_succcess_step[length(p_treatment_succcess_step)] else p_treatment_succcess_step[step + 1]
update(risk_multiplier_lar_nsp) = if (as.integer(step) >= length(risk_multiplier_lar_nsp_step)) risk_multiplier_lar_nsp_step[length(risk_multiplier_lar_nsp_step)] else risk_multiplier_lar_nsp_step[step+1]
update(gini_drug_deaths) = if (as.integer(step) >= length(gini_drug_deaths_step)) gini_drug_deaths_step[length(gini_drug_deaths_step)] else gini_drug_deaths_step[step+1]
update(immigrated_hcvpos_currstep) = if (as.integer(step) >= length(immigrated_hcvpos_step)) immigrated_hcvpos_step[length(immigrated_hcvpos_step)] else immigrated_hcvpos_step[step+1]
update(genpop_hcvpos_currstep) = if (as.integer(step) >= length(genpop_hcvpos_step)) genpop_hcvpos_step[length(genpop_hcvpos_step)] else genpop_hcvpos_step[step+1]
update(rate_treatment_input_currstep) = if (as.integer(step) >= length(rate_treatment_input_step)) rate_treatment_input_step[length(rate_treatment_input_step)] else rate_treatment_input_step[step+1]
update(rate_treatment_active_input_currstep) = if (as.integer(step) >= length(rate_treatment_active_input_step)) rate_treatment_active_input_step[length(rate_treatment_active_input_step)] else rate_treatment_active_input_step[step+1]
update(rate_death_active) = if (as.integer(step) >= length(rate_death_active_step)) rate_death_active_step[length(rate_death_active_step)] else rate_death_active_step[step+1]
update(rate_death_ex) = if (as.integer(step) >= length(rate_death_ex_step)) rate_death_ex_step[length(rate_death_ex_step)] else rate_death_ex_step[step+1]

# Update the time counter
update(time_check) = (step + 1) * dt 

# Outputs that are sanity checks
update(dt_check) = dt_check + 0 # Just for sanity checking
update(N_alive_PWID) = AN + AA + AC + AR + EN + EA + EC + ER + QN + QA + QC + QR
update(N_active) = AA + AN + AC + AR
update(N_ex_wr) = EA + EN + EC + ER
update(N_ex_wnr) = QA + QN + QC + QR
update(AA_ini_check) = AA_ini + 0 # Just for sanity checking, to get this variable as an output from the model
update(rate_treatment_rw_ini_check) = rate_treatment_rw_ini + 0 # Just for sanity checking, to get this variable as an output from the model
update(rate_debut_ini_check) = rate_debut_ini + 0 # Just for sanity checking, to get this variable as an output from the model
update(rate_IC_AC_check) = rate_IC_AC + 0 # Just for sanity checking, to get this variable as an output from the model

## Initial states:
initial(AN) = AN_ini - AA_ini - AC_ini # Removing those who are active and infected
initial(AA) = AA_ini
initial(AC) = AC_ini
initial(AR) = AR_ini
initial(EN) = EN_ini
initial(EA) = EA_ini
initial(EC) = EC_ini
initial(ER) = ER_ini
initial(QN) = QN_ini
initial(QA) = QA_ini
initial(QC) = QC_ini
initial(QR) = QR_ini
initial(D) = 0 # Cumulative deaths
initial(IC) = IC_ini
initial(IR) = IR_ini
initial(ID) = ID_ini
initial(GC) = GC_ini
initial(GR) = GR_ini
initial(GD) = GD_ini
# Rates:
initial(rate_debut) = rate_debut_ini 
initial(rate_quitting_wr) = rate_quitting_wr_ini
initial(rate_quitting_wnr) = rate_quitting_wnr_ini
initial(rate_relapsing) = rate_relapsing_ini
initial(beta) = beta_ini
initial(rate_recovery_ex) = rate_recovery_ex_ini
initial(rate_treatment) = rate_treatment_input_step[1]
initial(p_treatment_success) = p_treatment_succcess_step[1] 
initial(rate_treatment_active) = rate_treatment_active_input_step[1]
# initial(nsp_coverage) = nsp_coverage_step[1]
initial(risk_multiplier_lar_nsp) = risk_multiplier_lar_nsp_step[1]
initial(gini_drug_deaths) = gini_drug_deaths_step[1]
initial(immigrated_hcvpos_currstep) = immigrated_hcvpos_step[1]
initial(genpop_hcvpos_currstep) = genpop_hcvpos_step[1]
initial(rate_treatment_input_currstep) = rate_treatment_input_step[1]
initial(rate_treatment_active_input_currstep) = rate_treatment_active_input_step[1]
initial(rate_death_active) = rate_death_active_step[1]
initial(rate_death_ex) = rate_death_ex_step[1]
initial(switched_to_rw_treatment) = 0 # 0 = no, counting upward from 1 for each step from time_start_rw_treatment onwards
# Counters:
initial(cum_infections) = 0
initial(inc_infections) = 0
initial(cum_treatments) = 0
initial(inc_treatments_pwid) = 0
initial(inc_treatments_immigrants) = 0
initial(inc_treatments_genpop) = 0
initial(cum_debut) = 0
initial(inc_debut) = 0
initial(inc_death_active) = 0
initial(inc_death_ex_wr) = 0
initial(inc_death_ex_wnr) = 0
initial(inc_death_immigrants) = 0
initial(inc_death_genpop) = 0
initial(inc_immigrated_hcvpos) = 0
initial(cum_immigrated_hcvpos) = IC_ini
initial(cum_IC_AC) = 0
initial(cum_AA_AC) = 0
initial(time_check) = 0
initial(dt_check) = dt # Just for sanity checking
initial(N_active) = AA + AN + AC + AR # total active PWID population
initial(N_ex_wr) = EA + EN + EC + ER # total ex-PWID population. 
initial(N_ex_wnr) = QA + QN + QC + QR # total ex-PWID population. 
initial(N_alive_PWID) = AN + AA + AC + AR + EA + EN + EC + ER + QN + QA + QC + QR # total population.
initial(AA_ini_check) = AA_ini # Just for sanity checking
initial(rate_treatment_rw_ini_check) = rate_treatment_rw_ini # Just for sanity checking
initial(rate_debut_ini_check) = rate_debut_ini # Just for sanity checking
initial(rate_IC_AC_check) = rate_IC_AC # Just for sanity checking


## User defined parameters - default in parentheses:
rate_debut_ini = user(1)
rate_quitting_wr_ini = user(0.3)
rate_quitting_wnr_ini = user(0.3)
rate_relapsing_ini = user(0.3)
beta_ini = user(0.45)
sd_rw_debut = user(50) # Standard deviation on normal dist for random walk on PWID rate parameters
AN_ini = user(378) # From pwid data
AA_ini = user(0) # Seeding epidemic 
AC_ini = user(0) 
AR_ini = user(0)
EA_ini = user(0) 
EN_ini = user(35) # From pwid data
EC_ini = user(0)
ER_ini = user(0) 
QA_ini = user(0) 
QN_ini = user(0) # From pwid data
QC_ini = user(0)
QR_ini = user(0) 
IC_ini = user(0)
IR_ini = user(0)
ID_ini = user(0)
GC_ini = user(0)
GR_ini = user(0)
GD_ini = user(0)
rate_recovery_ex_ini = user(0.26*(1/0.5)) # Recovery rate due to natural clearance of infection. Assuming there's a 26 percent probability of clearing infection during 6 months.
rel_prob_recovery_active = user(1.0) # Parameter to modulate recovery probability of active PWID compared to ex-PWID
rel_prob_reinfection = user(1.0) # Parameter to modulate reinfection probability 
dt = user(0.25) # Timestep
rate_acute_chronic = user(0.74*1/(0.5)) # (2 yr^-1) -- 6 months average time for infection to become chronic (0.5 yr^-1) * 74% probability of becoming chronic
rate_treatment_rw_ini = user(1)
rate_treatment_input_step[] = user() # Rate of treatment initiation, for the period when treatments should not be random walked
dim(rate_treatment_input_step) = user()
rate_treatment_active_input_step[] = user()
dim(rate_treatment_active_input_step) = user()
time_start_rw_treatment = user() # Time when random walk on treatment rate should start
sd_rtreat = user(1) 
rate_treatment_max = user(0.8) # Maximum rate of treatment initiation, used to set upper limit on random walk
rate_death_active_step[] = user() # Dominated by overdose deaths. Taken to be time-varying input
dim(rate_death_active_step) = user()
rate_death_ex_step[] = user() # Natural age-related death rate for ex-PWID
dim(rate_death_ex_step) = user()
risk_multiplier_lar_nsp_step[] = user() # Risk multiplier for PWID using LAR and/or NSP, modulating transmission rate
dim(risk_multiplier_lar_nsp_step) = user()
gini_drug_deaths_step[] = user()
dim(gini_drug_deaths_step) = user()
p_treatment_succcess_step[] = user() # Vector of year-varying probability of successful treatment
dim(p_treatment_succcess_step) = user()
# Immigrants and general population:
rate_IC_AC = user(0)
rate_death_imm = user(0)
immigrated_hcvpos_step[] = user()
dim(immigrated_hcvpos_step) = user()
rate_death_genpop = user(0)
genpop_hcvpos_step[] = user()
dim(genpop_hcvpos_step) = user()
# Other
time_forecast = user() # Time when forecast starts, when the model should be run with no random walks
