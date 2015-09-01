function wind = GenWind(f_mean,f_sigma,l_mean,l_sigma,total_time,dt)

f_series = normrnd(f_mean,f_sigma,50000,1); 
l_series = normrnd(l_mean,l_sigma,50000,1);

l_series(l_series<0)=0;

if sum(l_series)<total_time
        error('adapth time parameters')
end


i=1;
k=1;
wind=zeros(round(sum(l_series/dt)),1);
while i<(total_time/dt)
    step_nrs = round(l_series(k)/dt);
    wind(i:(i+step_nrs))=f_series(k);
    k=k+1;
    i=i+step_nrs;
end
%crop wind
wind = wind(1:total_time/dt,1);
wind(wind<0)=0;

    









