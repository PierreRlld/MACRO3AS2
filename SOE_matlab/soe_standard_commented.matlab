%----------------------------------------------------------------
close all;
%format long
%----------------------------------------------------------------

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var c_H r_H pic_H pi_H mc_H w_H h_H y_H p_H NFA_H lb_H ex_H
	c_F r_F pic_F pi_F mc_F w_F h_F y_F p_F NFA_F lb_F ex_F
	de rer
	e_H e_F tau_H tau_F mu_H mu_F g_H g_F varrho_H varrho_F
	e_z_H e_p_H e_r_H e_x_H e_t_H e_g_H
	e_z_F e_p_F e_r_F e_x_F e_t_F e_g_F 
	gy_H_obs gy_F_obs gc_H_obs gc_F_obs pi_H_obs pi_F_obs r_H_obs r_F_obs de_obs  drer_obs ex_H_obs ex_F_obs
	e_e;

varexo	eta_z_H eta_p_H eta_r_H eta_x_H eta_z_F eta_p_F eta_r_F eta_x_F eta_e eta_t_H eta_t_F  eta_g_H eta_g_F ;

parameters	sigmaC_H sigmaC_F sigmaH_H sigmaH_F beta alpha hc_H hc_F chi_B chi_H chi_F xi_H xi_F epsilon_H epsilon_F mu alphaC_H alphaC_F n rho phi_pi phi_y 
			rho_e phi_H phi_F tau0_H tau0_F y0 sig_H sig_F theta1 theta2 varphi A psi piss Hss gy_H gy_F
			rho_z_H rho_r_H rho_p_H rho_x_H rho_t_H rho_g_H
			rho_z_F rho_r_F rho_p_F rho_x_F rho_t_F rho_g_F
			;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

% - - - - - - - - - - - - - - -
%% Difference between H and F :
%% curvatures sigmaC and sigmaH
%% consumption habits hc_H and hc_F
%% cost of adjusting price is different for H and F firms
%% substitution between goods for retailers is different
%% 
% - - - - - - - - - - - - - - -

sigmaC_H		= 1.5;		% risk aversion
sigmaC_F		= 1.2;		% risk aversion
sigmaH_H		= 2;		% labor supply
sigmaH_F		= 1.9;		% labor supply
beta		= .994;			% discount factor
alpha		= .7;			% share of labor in production
hc_H		= .7;			% consumption habits
hc_F		= .6;			% consumption habits
chi_B		= 0.007;		% cost of foreign debt (same for hh in H or in F)
xi_H		= 100;			% cost adjustment prices
xi_F		= 80;			% cost adjustment prices
epsilon_H	= 10;			% imperfect substitution between goods
epsilon_F	= 9.5;			% imperfect substitution between goods
mu			= 2;			% Substitution between home/foreign goods
alphaC_H	= .1;			% Share of home goods in consumption basket
alphaC_F	= .11;			% Share of foreign? goods in consumption basket
rho			= .8;			% Monetary policy coefficient smoothing
phi_pi		= 1.5;			% Monetary policy reaction to inflation
phi_y		= .05;			% Monetary policy reaction to output
n			= .4;			% share of home country then size of foreign country 1-n
varphi		= 0.2;			% elasticity of emission to GDP
piss		= 1.005;		% steady state inflation
gy_H 		= 0.2;			% Public spending to gdp
gy_F 		= 0.2;			% Public spending to gdp

% value of main variables:
tau0_H	= 50 /1000;	% value of carbon tax ($/ton)
tau0_F	= 50 /1000;	% value of carbon tax ($/ton)
sig_H	= 0.2; 		% Carbon intensity USA 0.2 Gt / Trillions USD
sig_F	= 0.2; 		% Carbon intensity USA 0.2 Gt / Trillions USD
y0	 	= 25;		% trillions usd PPA https://data.worldbank.org/indicator/NY.GDP.MKTP.CD
theta1  = 0.05;		% level of abatement costs
theta2  = 2.6;		% curvature abatement cost
Hss		= 1/3;		% labor supply in ss


% parameters of autoregressive shocks						
rho_z_H 	= .95; rho_z_F 	= .95;
rho_p_H 	= .95; rho_p_F 	= .95;
rho_r_H		= .4;  rho_r_F	= .4;
rho_x_H		= .4;  rho_x_F	= .4;
rho_t_H		= .4;  rho_t_F	= .4;
rho_g_H		= .4;  rho_g_F	= .8;
rho_e		= .1;


%% SS : CLOSE-FORM SOLUTION OF THE STEADY-STATE
steady_state_model;
	h_H		= Hss;
	h_F		= Hss;
	tau_H 	= tau0_H;
	tau_F 	= tau0_F;
	phi_H	= (1-(1-n)*alphaC_H);
	phi_F	= (1-n*alphaC_F);
	y_H     = y0/n;
	A		= y_H/h_H^alpha;
	y_F		= A*h_F^alpha;
	g_H     = gy_H*y_H;
	g_F     = gy_H*y_F;
	mu_H	= (tau_H*sig_H*y_H^(1-varphi)/(theta2*theta1))^(1/(theta2-1));
	mu_F	= (tau_F*sig_F*y_F^(1-varphi)/(theta2*theta1))^(1/(theta2-1));
	e_H 	= n*sig_H*(1-mu_H)*y_H^(1-varphi);
	e_F 	= (1-n)*sig_H*(1-mu_F)*y_F^(1-varphi);
	c_F 	= (y_F - g_F - theta1*mu_F^theta2*y_F -((1-phi_H)*n/(1-n))* (y_H-g_H-theta1*mu_H^theta2*y_H)/phi_H)/((phi_F-(1-phi_H)*(1-phi_F)/phi_H));
	c_H		= (y_H-g_H-theta1*mu_H^theta2*y_H-(1-phi_F)*c_F*(1-n)/n)/phi_H  ;
	lb_H 	= (c_H-hc_H*c_H)^-sigmaC_H;
	lb_F 	= (c_F-hc_F*c_F)^-sigmaC_F;
	r_H		= piss/beta;
	r_F		= piss/beta;
	pic_H	= piss; pic_F	= piss;
	pi_H	= piss; pi_F	= piss;
	de		= 1;
	NFA_H	= ((phi_H*c_H + (1-phi_F)*c_F*(1-n)/n) - c_H)/(1-r_F/pic_H);
	NFA_F	= -n/(1-n)*NFA_H;
	rer		= 1;
	mc_H	= (epsilon_H-1)/epsilon_H;
	mc_F	= (epsilon_F-1)/epsilon_F;
	varrho_H = mc_H - theta1*mu_H^theta2 - tau_H*(1-varphi)*sig_H*(1-mu_H)*y_H^(-varphi);
	varrho_F = mc_F - theta1*mu_F^theta2 - tau_F*(1-varphi)*sig_F*(1-mu_F)*y_F^(-varphi);
	w_H		= varrho_H/h_H*(alpha*y_H);
	w_F		= varrho_F/h_F*(alpha*y_F);
	p_H		= 1; 
	p_F		= 1;
	ex_H 	= (1-phi_F)*c_F*(1-n);
	ex_F 	= (1-phi_H)*c_H*n;
	chi_H	= lb_H*w_H/(h_H^sigmaH_H);
	chi_F	= lb_F*w_F/(h_F^sigmaH_F);
	e_z_H	= 1; e_p_H = 1; e_r_H = 1; e_x_H = 1; e_t_H = 1; e_g_H  = 1; e_e = 1;
	e_z_F	= 1; e_p_F = 1; e_r_F = 1; e_x_F = 1; e_t_F = 1; e_g_F  = 1;
	gy_H_obs = 0; gy_F_obs = 0; gc_H_obs = 0; gc_F_obs = 0; pi_H_obs = 0; pi_F_obs = 0; r_H_obs = 0; r_F_obs = 0; de_obs = 0;  drer_obs = 0; ex_H_obs = 0; ex_F_obs = 0;
end;

%----------------------------------------------------------------
% 3. Model (the number refers to the equation in the paper)
%----------------------------------------------------------------

model;
	% ==============
	%%% Households
	[name='FOC c']
	lb_H = (c_H-hc_H*c_H(-1))^-sigmaC_H;
	lb_F = (c_F-hc_F*c_F(-1))^-sigmaC_F;

	[name='Euler equation']
	lb_H = beta*lb_H(+1)*r_H/pic_H(+1);
	lb_F = beta*lb_F(+1)*r_F/pic_F(+1);

	[name='Labor Supply']
	chi_H*h_H^sigmaH_H = lb_H*w_H;
	chi_F*h_F^sigmaH_F = lb_F*w_F;
	
	%%% FIRMS
	[name='NKPC']
	(1-epsilon_H) + epsilon_H*e_p_H*mc_H - xi_H*pi_H*(pi_H-piss) + xi_H*beta*((c_H(+1)-hc_H*c_H)/(c_H-hc_H*c_H(-1)))^-sigmaC_H*pi_H(+1)*(pi_H(+1)-piss)*y_H(+1)/y_H; % (sous entendu = 0)
 	(1-epsilon_F) + epsilon_F*e_p_F*mc_F - xi_F*pi_F*(pi_F-piss) + xi_F*beta*((c_F(+1)-hc_F*c_F)/(c_F-hc_F*c_F(-1)))^-sigmaC_F*pi_F(+1)*(pi_F(+1)-piss)*y_F(+1)/y_H; % (sous entendu = 0)

	[name='FOC h']
	varrho_H = h_H/(alpha*y_H)*w_H;
	varrho_F = h_F/(alpha*y_F)*w_F;

	[name='Production function']
	y_H = A*e_z_H*h_H^alpha;
	y_F = A*e_z_F*h_F^alpha;

	[name='CES price index']
	1 = phi_H*p_H^(1-mu) + (1-phi_H)*rer^(1-mu);
	1 = phi_F*p_F^(1-mu) + (1-phi_F)*(1/rer)^(1-mu);

	[name='Relative price']
	p_H/p_H(-1) = pi_H/pic_H;
	p_F/p_F(-1) = pi_F/pic_F;

	[name='Total emissions']
	e_H = n*sig_H*(1-mu_H)*y_H^(1-varphi);
	e_F = (1-n)*sig_F*(1-mu_F)*y_F^(1-varphi);

	[name='FOC mu']
	varrho_H = mc_H - theta1*mu_H^theta2 - tau_H*(1-varphi)*sig_H*(1-mu_H)*y_H^-varphi;
	varrho_F = mc_F - theta1*mu_F^theta2 - tau_F*(1-varphi)*sig_F*(1-mu_F)*y_F^-varphi;

	[name='FOC y']
	tau_H*sig_H*y_H^(1-varphi) = theta2*theta1*mu_H^(theta2-1);
	tau_F*sig_F*y_F^(1-varphi) = theta2*theta1*mu_F^(theta2-1);
	% ==============


	% ==============
	%%% AGGREGATION
	[name='Resources constraint']
	y_H = phi_H*p_H^-mu*c_H + e_x_F*(1-phi_F)*(p_H/rer)^-mu*c_F*(1-n)/n  + g_H + theta1*mu_H^theta2*y_H + 0.5*xi_H*(pi_H-piss)^2*y_H + 0.5*chi_B*(NFA_H-STEADY_STATE(NFA_H))^2;
	y_F = phi_F*p_F^-mu*c_F + e_x_H*(1-phi_H)*(p_F*rer)^-mu*c_H*n/(1-n)  + g_F + theta1*mu_F^theta2*y_F + 0.5*xi_F*(pi_F-piss)^2*y_F - 0.5*chi_B*(NFA_F-STEADY_STATE(NFA_F))^2;
	% ==============
	

	% ==============
	%%% POLICIES
	[name='Monetary Policy Rule']
	r_H = r_H(-1)^rho * (STEADY_STATE(r_H)*(pic_H/steady_STATE(pic_H))^phi_pi*(y_H/STEADY_STATE(y_H))^phi_y)^(1-rho)*e_r_H;
	r_F = r_F(-1)^rho * (STEADY_STATE(r_F)*(pic_F/steady_STATE(pic_F))^phi_pi*(y_F/STEADY_STATE(y_F))^phi_y)^(1-rho)*e_r_F;

	[name='Public spending']
	g_H = gy_H*steady_state(y_H)*e_g_H;
	g_F = gy_F*steady_state(y_F)*e_g_F;

	[name='Carbon tax']
	tau_H = tau0_H*e_t_H;
	tau_F = tau0_F*e_t_F;
	% ==============


	% ==============
	%%% Common macro variables from the Home country perspective
	[name='Net Foreign assets accumulation']
	NFA_H = r_F(-1)/pic_H*NFA_H(-1)*de + p_H*(phi_H*p_H^-mu*c_H + e_x_F*(1-phi_F)*(1-n)/n*(p_H/rer)^-mu*c_F) - c_H;

	[name='International financial markets accounting']
	n*NFA_H + (1-n)*NFA_F = 0;

	[name='Nominal exchange rate growth']
	de(+1) = (1+chi_B*(NFA_H-STEADY_STATE(NFA_H)))*r_H/r_F/e_e;

	[name='Real exchange rate']
	rer/rer(-1) = de*pi_F/pi_H;

	[name='Exports']
	ex_H = e_x_H*(1-phi_F)*(p_H/rer)^-mu*c_F*(1-n);
	ex_F = e_x_F*(1-phi_H)*(p_F*rer)^-mu*c_H*n;
	% ==============
	

	% ==============
	%% Observable variables > data from dbnomics should have EXACT SAME names for the chosen series
	[name='measurement GDP']
	gy_H_obs = log(y_H/y_H(-1));
	gy_F_obs = log(y_F/y_F(-1));

	[name='measurement consumption']
	gc_H_obs = log(c_H/c_H(-1));
	gc_F_obs = log(c_H/c_H(-1));

	[name='measurement inflation']
	pi_H_obs = pi_H - steady_state(pi_H);
	pi_F_obs = pi_F - steady_state(pi_F);

	[name='measurement interest rate']
	r_H_obs  = r_H  - steady_state(r_H);
	r_F_obs  = r_F  - steady_state(r_F);

	[name='measurement nominal exchange rate change']
	de_obs  = log(de);

	[name='measurement real exchange rate change']
	drer_obs  = log(rer/rer(-1));

	[name='measurement exports change']
	ex_H_obs  = log(ex_H/ex_H(-1));
	ex_F_obs  = log(ex_F/ex_F(-1));
	% ==============


	% ==============
	%% Stochastic processes
	[name='Country specific shocks']
	%% Home
	log(e_z_H) = rho_z_H*log(e_z_H(-1)) + eta_z_H;	
	log(e_p_H) = rho_p_H*log(e_p_H(-1)) + eta_p_H;
	log(e_r_H) = rho_r_H*log(e_r_H(-1)) + eta_r_H;
	log(e_x_H) = rho_x_H*log(e_x_H(-1)) + eta_x_H;
	log(e_g_H) = rho_x_H*log(e_g_H(-1)) + eta_g_H;
	log(e_t_H) = rho_t_H*log(e_t_H(-1)) + eta_t_H;
	%% Foreign
	log(e_z_F) = rho_z_F*log(e_z_F(-1)) + eta_z_F;
	log(e_p_F) = rho_p_F*log(e_p_F(-1)) + eta_p_F;
	log(e_r_F) = rho_r_F*log(e_r_F(-1)) + eta_r_F;
	log(e_x_F) = rho_x_F*log(e_x_F(-1)) + eta_x_F;
	log(e_g_F) = rho_g_F*log(e_g_F(-1)) + eta_g_F;
	log(e_t_F) = rho_t_F*log(e_t_F(-1)) + eta_t_F;

	log(e_e)   = rho_e*log(e_e(-1)) + eta_e;
	
end;



% > check the starting values for the steady state
resid;

% > compute steady state given the starting values
steady;

% > check Blanchard-Kahn-conditions :
check;


%%% Stochastic Simulations // replace with your codes
shocks;
var eta_z_H;  stderr 0.01;
var eta_p_H;  stderr 0.01;
var eta_r_H;  stderr 0.01;
var eta_e;	  stderr 0.01;
var eta_x_H;  stderr 0.01;
end;

% Basic simulation of the model :
stoch_simul(order=1, irf=20) y_H y_F c_H c_F pi_H pi_F r_H r_F rer ex_H ex_F;


%----------------------------------------------------------------
% 4. Estimation
%----------------------------------------------------------------
%% Now we have to : estimate parameters to mimic a few observed variables > measurement equations
%% evaluate the 'estimated' model

% ----------------------------
% (1) Provide observable series :
% Dynare will look for these exact variables in the data matrix
% varobs gy_obs pi_obs r_obs gc_obs gi_obs;
% ----------------------------
