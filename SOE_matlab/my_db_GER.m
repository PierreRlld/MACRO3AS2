[output_mat,output_table,dates_nb] =  call_dbnomics('OECD/QNA/DEU.B1_GE.CQRSA.Q','OECD/QNA/DEU.B1_GE.DNBSA.Q','IMF/DOT/Q.DE.TMG_CIF_USD.US','ECB/EXR/Q.USD.EUR.SP00.A','ECB/EXR/Q.E01.USD.ERC0.A');

%% ================================
% --- Observed series ---
% (1) [GDP€] : "Germany – Gross domestic product - expenditure approach – National currency, current prices, quarterly levels, seasonally adjusted – Quarterly" ['OECD/QNA/DEU.B1_GE.CQRSA.Q']
% (1.1) [GDP Deflator 2015 = 100] : "Germany – Gross domestic product - expenditure approach – Deflator, national base/reference year, seasonally adjusted – Quarterly" ['OECD/QNA/DEU.B1_GE.DNBSA.Q']
% (2) [Imports from US in $ = US exports to GER] : "Quarterly – Germany – Goods, Value of Imports, Cost, Insurance, Freight (CIF), US Dollars – United States, Millions" ['IMF/DOT/Q.DE.TMG_CIF_USD.US']
% (3) [Nominal exchange rate 1€=x$] : "Quarterly – US dollar – Euro – Spot – Average" [ECB/EXR/Q.USD.EUR.SP00.A]
% (4) 
% (?) [Real Exchange Rate (Index)] : "Quarterly – Narrow EER group of trading partners (fixed composition) – US dollar – Real effective exch. rate CPI deflated – Average" ['ECB/EXR/Q.E01.USD.ERC0.A']


%% ================================
idx = find(~isnan(sum(output_mat(:,2:end),2)));	% Find rows without NaN
data = output_mat(idx,:);						% > dataset without missing values
T = dates_nb(idx);								% > corresponding dates
imports_eur = data(:,4)./data(:,5);             % Imports are in $ : convert in € using NER

% cpi quarterly ???????????

%% ================================
% --- Define "observed series" ---
% (1) Real output growth rate
gy_H_obs = diff(log(output_mat(:,1)));
% > Measurement equation : gy_H_obs = log(y_H/y_H(-1))

% (2) Foreign exports = Home imports
ex_F_obs = diff(log(output_mat(:,2)));
% > Measurement equation : ex_F_obs  = log(ex_F/ex_F(-1))

% (3) Nominal exchange rate change (growth rate) : Foreign->Home ($->€)
% drer_obs = diff(log(output_mat(:,3)));
% > Measurement equation : de_obs = log(de) (où de = ratio des exchange rate dans modèle)


%% ================================
% --- Export ---
% save into `ger_obs` the series selected observed series
% save ger_obs ;



