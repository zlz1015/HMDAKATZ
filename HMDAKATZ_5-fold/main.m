clear;
for cv=1:10
    %getParameters();
    HMDAKATZ_5CV(1,1,0.1,2)
    overallauc(cv)=positiontooverallauc();
end
save overallauc overallauc
a=mean(overallauc)
b=std(overallauc)

