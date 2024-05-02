close all;

%% >>> Load estimated parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("soe_res_save.mat");

% load ESTIMATED PARAMETERS
fn = fieldnames(oo_.posterior_mean.parameters);
for ix = 1:size(fn,1)
	set_param_value(fn{ix},eval(['oo_.posterior_mean.parameters.' fn{ix} ]))
end
% load ESTIMATED SHOCKS
fx = fieldnames(oo_.posterior_mean.shocks_std);
for ix = 1:size(fx,1)
	idx = strmatch(fx{ix},M_.exo_names,'exact');
	M_.Sigma_e(idx,idx) = eval(['oo_.posterior_mean.shocks_std.' fx{ix}])^2;
end

load(options_.datafile);
% dataset_ object stores obs series
if exist('T') ==1
	Tvec = T;
else
	Tvec = 1:size(dataset_,1);
end
Tfreq = mean(diff(Tvec));


%% >>> END OF SAMPLE FORECASTING - PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tprior = 20; 	% period before forecasts to plot
Tvec2 = Tvec(end) + (0:(options_.forecast))*Tfreq;
for i1 = 1 :size(dataset_.name,1)
	% position indices of observed variable number i1 in dataset and Model 
	idv		= strmatch(dataset_.name{i1},M_.endo_names,'exact');
	idd		= strmatch(dataset_.name{i1},dataset_.name,'exact');
	if ~isempty(idd) && isfield(oo_.MeanForecast.Mean, dataset_.name{i1})
		% Find model based values for observed variables (SmoothedVariables) + add observed mean
		% because the model is based on zero-mean version of the data (prefilter in estimation)
		yobs   = eval(['oo_.SmoothedVariables.' dataset_.name{i1}])+dataset_info.descriptive.mean(idd);
		yfc    = eval(['oo_.MeanForecast.Mean.'  dataset_.name{i1}])+dataset_info.descriptive.mean(idd);
		yfcVar = sqrt(eval(['oo_.MeanForecast.Var.' dataset_.name{i1}]));
        % >> Create individual chart :
		figure;
		plot(Tvec(end-tprior+1:end),yobs(end-tprior+1:end))
		hold on;
			plot(Tvec2,[yobs(end) yfc'] ,'r--','LineWidth',1.5);
			% '1.96' due to normality assumption for 'residuals' in space-state representation of the model
			plot(Tvec2,[yobs(end) (yfc+1.96*yfcVar)'],'r:','LineWidth',1.5)
			plot(Tvec2,[yobs(end) (yfc-1.96*yfcVar)'],'r:','LineWidth',1.5)
			grid on;
			xlim([Tvec(end-tprior+1) Tvec2(end)])
			legend('Sample','Forecasting','Uncertainty')
			title(['forecasting of ' M_.endo_names_tex{idv}])
		hold off;
	else
		warning([ dataset_.name{i1} ' is not an observable or you have not computed its forecast'])
	end
end


%% >>> COUNTERFACTUAL EXERCISES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stack estimated values for exogenous shocks in a matrix
fx = fieldnames(oo_.SmoothedShocks);
for ix=1:size(fx,1)
	% extract the correct (model-based) series from oo_.SmoothedShocks
	shock_mat = eval(['oo_.SmoothedShocks.' fx{ix}]);
	if ix==1; ee_mat = zeros(length(shock_mat),M_.exo_nbr); end;
	ee_mat(:,strmatch(fx{ix},M_.exo_names,'exact')) = shock_mat;
end

% ------
%>>> Simulate BASELINE scenario
% SOLVE DECISION RULEs
[oo_.dr, info, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);
% SIMULATE the model
y_            = simult_(M_,options_,oo_.dr.ys,oo_.dr,ee_mat,options_.order);

% ------
%>>> Simulate ALTERNATIVE scenario
% (!) make a copy (!)
Mx  = M_;
oox = oo_;
% (!) CHANGE PARAMETER (!)
Mx.params(strcmp('phi_y',M_.param_names)) = .25;
% solve new decision rule
[oox.dr, info, Mx.params] = resol(0, Mx, options_, oox.dr, oox.dr.ys, oox.exo_steady_state, oox.exo_det_steady_state);
% simulate dovish central bank
ydov            = simult_(Mx,options_,oox.dr.ys,oox.dr,ee_mat,options_.order);


% ------
% Plot results
var_names={'lny','lnc','lni','lnpi','lnr','h_obs'};
Ty = [T(1)-Tfreq;T];
% draw_tables from 'draw_tables.m'
draw_tables(var_names,M_,Ty,[],y_,ydov)
legend('Estimated','Dovish')