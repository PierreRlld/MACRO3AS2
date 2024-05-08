%%%%%%%%%%%%%%%%% END OF SAMPLE FORECASTING - PLOTS
tprior = 20; % period before forecasts to plot
Tvec2 = Tvec(end) + (0:(options_.forecast))*Tfreq;
for i1 = 1 :size(dataset_.name,1)
	idv		= strmatch(dataset_.name{i1},M_.endo_names,'exact');
	idd		= strmatch(dataset_.name{i1},dataset_.name,'exact');
	if ~isempty(idd) && isfield(oo_.MeanForecast.Mean, dataset_.name{i1})
		% Create chart :
		% Find model based values for observed variables (SmoothedVariables) + add observed mean
		% because the model is based on zero-mean version of the data (prefilter in estimation)
		yobs   = eval(['oo_.SmoothedVariables.' dataset_.name{i1}])+dataset_info.descriptive.mean(idd);
		yfc    = eval(['oo_.MeanForecast.Mean.'  dataset_.name{i1}])+dataset_info.descriptive.mean(idd);
		yfcVar = sqrt(eval(['oo_.MeanForecast.Var.' dataset_.name{i1}]));
		figure;
		plot(Tvec(end-tprior+1:end),yobs(end-tprior+1:end))
		hold on;
			plot(Tvec2,[yobs(end) yfc'] ,'r--','LineWidth',1.5);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%% COUNTERFACTUAL EXERCISES %%%%%%%%%%%%%%%%%%
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
% make a copy
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
var_names={'lny_H','lnc_H','lnpi_H','lnr_H'};
Ty = [T(1)-Tfreq;T];
draw_tables(var_names,M_,Ty,[],y_,ydov)
legend('Estimated','Dovish')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%% FORECAST UNDER ALTERNATIVE POLICY %%%%%%%%%%%%%%%%%%
Thorizon 	= 12; % number of quarters for simulation
% Built baseline forecast
fx = fieldnames(oo_.SmoothedShocks);
for ix=1:size(fx,1)
	shock_mat = eval(['oo_.SmoothedShocks.' fx{ix}]);
	if ix==1; ee_mat2 = zeros(length(shock_mat),M_.exo_nbr); end;
	ee_mat2(:,strmatch(fx{ix},M_.exo_names,'exact')) = shock_mat;
end
% add mean-wise forecast with zero mean shocks
ee_mat2 	= [ee_mat;zeros(Thorizon,M_.exo_nbr)];
Tvec2 		= Tvec(1):Tfreq:(Tvec(1)+size(ee_mat2,1)*Tfreq);

% ------
%>>> Simulate ALTERNATIVE scenario
% SOLVE DECISION RULEs
[oo_.dr, info, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);
% SIMULATE the model
y_            = simult_(M_,options_,oo_.dr.ys,oo_.dr,ee_mat2,options_.order);

%%% Add a positive fiscal shock
% make a copy of shock matrix
ee_matx = ee_mat2;
% select fiscal shock
idx = strmatch('eta_g_H',M_.exo_names,'exact');
ee_matx(end-Thorizon+1,idx) = 0.05;% add a 5 percent increase in public spending
% simulate the model
y_fiscal           = simult_(M_,options_,oo_.dr.ys,oo_.dr,ee_matx,options_.order);


%%% Add a positive carbon shock
% make a copy of shock matrix
ee_matx = ee_mat2;
% select fiscal shock
idx = strmatch('eta_t_H',M_.exo_names,'exact');
ee_matx(end-Thorizon+1,idx) = 0.5;% add a 50 percent increase in carbon price 
% simulate the model
y_carbon           = simult_(M_,options_,oo_.dr.ys,oo_.dr,ee_matx,options_.order);

%%% Add a negative monetary policy shock
% make a copy of shock matrix
ee_matx = ee_mat2;
% select fiscal shock
idx = strmatch('eta_r_H',M_.exo_names,'exact');
ee_matx(end-Thorizon+1,idx) = -0.05;% add a 50 percent increase in carbon price 
% simulate the model
y_monetary           = simult_(M_,options_,oo_.dr.ys,oo_.dr,ee_matx,options_.order);

% draw result
var_names={'lny_H','lnpi_H','lnr_H','g_H','tau_H'};
Ty = [T(1)-Tfreq;T];
draw_tables(var_names,M_,Tvec2,[2023 Tvec2(end)],y_,y_fiscal,y_carbon,y_monetary)
legend('Estimated','Fiscal','Carbon','Monetary')