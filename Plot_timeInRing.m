black.col = {'k'};
blue.col = {'b'};
green.col = {'g'};
red.col = {'r'};
mseb(radii,nanmean((R(:,prevT0==0 & bait0==0)),2)',nanstderr((R(:,prevT0==0 & bait0==0)),[],2)',black);
mseb(radii,nanmean((R(:,prevT0==1 & bait0==0)),2)',nanstderr((R(:,prevT0==1 & bait0==0)),[],2)',blue);
mseb(radii,nanmean((R(:,prevT0==0 & bait0==1)),2)',nanstderr((R(:,prevT0==0 & bait0==1)),[],2)',green);
mseb(radii,nanmean((R(:,prevT0==1 & bait0==1)),2)',nanstderr((R(:,prevT0==1 & bait0==1)),[],2)',red);
xlabel('Distance from Spot (cm)','fontsize',24);
ylabel('time in ring / area (s/cm^2)','fontsize',24);
legend('Prev = unbaited, Curr = unbaited','Prev = baited, Curr = unbaited','Prev = unbaited, Curr = baited','Prev = baited, Curr = baited');
