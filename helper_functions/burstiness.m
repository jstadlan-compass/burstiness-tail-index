function BP = burstiness(data)
    CoV = std(data)./mean(data);
    BP = (CoV-1)./(CoV+1);

