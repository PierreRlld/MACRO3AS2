[output_table,~,T] = call_dbnomics('OECD/QNA/DEU.B1_GE.CQRSA.Q','OECD/QNA/DEU.B1_GE.DNBSA.Q','OECD/KEI/IR3TIB01.DEU.ST.Q','OECD/QNA/DEU.P71.CQRSA.Q');

%% ================================
% --- Observed series ---
% (1) [GDP€] : "Germany – Gross domestic product - expenditure approach – National currency, current prices, quarterly levels, seasonally adjusted – Quarterly" ['OECD/QNA/DEU.B1_GE.CQRSA.Q']
% (2) [GDP Deflator 2015 = 100] : "Germany – Gross domestic product - expenditure approach – Deflator, national base/reference year, seasonally adjusted – Quarterly" ['OECD/QNA/DEU.B1_GE.DNBSA.Q']
% (3) [Nominal interest rate] : "3 month interbank rate – Germany" ['OECD/KEI/IR3TIB01.DEU.ST.Q']
% (4) [German imports = Foreign exports] : "Germany – Imports of goods – National currency, current prices, quarterly levels, seasonally adjusted – Quarterly" ['OECD/QNA/DEU.P71.CQRSA.Q']
		% goods only : ['OECD/QNA/DEU.P71.CQRSA.Q'] ?
		% goods and services : ['OECD/QNA/DEU.P7.CQRSA.Q'] ?

% remark 9 chained values 

%% ================================

% select non NaN ids
idx 			= find(~isnan(sum(output_table(:,2:end),2))); % Find rows without NaN
output_table 	= output_table(idx,:);                        % > dataset without missing values
T				= T(idx);                                     % > corresponding dates


% > Normalize prices to 1 in 2015
id2015 = find(T==2015);
def = output_table(:,3)/output_table(id2015,3);

%% ================================
% --- Define "observed series" ---
% Should demean the sample to make it consistant with model definition of
% deviations from steady-state. We do not do it here but rather in the
% 'estimation' part of soe.mod with 'prefilter=1'.

% (1) Real output growth :
gy_H_obs  = diff(log(output_table(:,2)./(def)));
% > Measurement equation : gy_H_obs = log(y_H/y_H(-1))
% > Where y_h = model total output measured in goods

% (2) Inflation :
pi_H_obs  = diff(log(def));
% > Measurement equation : log(p_H/p_H(-1))

% (3) Nominal interest rate :
% Rate data are annualized and *100
% Remove one observation to be consistant with other variables that lose
% one due to differentiation.
r_H_obs	= output_table(2:end,4)/400;

% (4) Real Foreign exports growth :
ex_F_obs  = diff(log(output_table(:,5)./(def)));
% > Measurement equation : ex_F_obs  = log(ex_F/ex_F(-1))
% > Where ex_F = model H imports of Foreign products measured in goods

T = T(2:end);

% save into myobsGER.mat
save myobsGER gy_H_obs T pi_H_obs r_H_obs ex_F_obs;

figure;

subplot(2,2,1)
plot(T,gy_H_obs)
xlim([min(T) max(T)]);
title('Real output (%QoQ)')

subplot(2,2,2)
plot(T,ex_F_obs)
xlim([min(T) max(T)]);
title('German imports (%QoQ)')

subplot(2,2,3)
plot(T,pi_H_obs)
xlim([min(T) max(T)]);
title('Inflation (%QoQ)')

subplot(2,2,4)
plot(T,r_H_obs)
xlim([min(T) max(T)]);
title('Nominal Interest rate')