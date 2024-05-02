[output_table,~,T] = call_dbnomics('OECD/QNA/DEU.B1_GE.CQRSA.Q','OECD/QNA/DEU.B1_GE.DNBSA.Q','OECD/KEI/IR3TIB01.DEU.ST.Q','OECD/QNA/DEU.P7.CQRSA.Q');
% Output (D), Deflator (D), Nominal Rate (D), Exports (F)

% select non NaN ids
idx 			= find(~isnan(sum(output_table(:,2:end),2)));
output_table 	= output_table(idx,:);
T				= T(idx);


% we normalize to one prices and in population for 2015
id2015 = find(T==2015);
def = output_table(:,3)/output_table(id2015,3);


%% taking in real growth rates per capita
gy_H_obs  = diff(log(output_table(:,2)./(def)));
ex_F_obs  = diff(log(output_table(:,5)./(def)));

% inflation rate
pi_H_obs  = diff(log(def));
% quarterly interest rate
r_H_obs	= output_table(2:end,4)/400;

T = T(2:end);

% save into myobsGER.mat
save myobsGER gy_H_obs T pi_H_obs r_H_obs ex_F_obs;

figure;
subplot(2,2,1)
plot(T,gy_H_obs)
xlim([min(T) max(T)]);
title('output growth')
subplot(2,2,2)
plot(T,ex_F_obs)
xlim([min(T) max(T)]);
title('imports growth')
subplot(2,2,3)
plot(T,pi_H_obs)
xlim([min(T) max(T)]);
title('inflation rate')
subplot(2,2,4)
plot(T,r_H_obs)
xlim([min(T) max(T)]);
title('nominal rate')