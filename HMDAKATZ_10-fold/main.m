clear;
for cv=1:100
    HMDAKATZ_10CV(1,1,0.1,2);
    overallauc(cv)=positiontooverallauc();
end
save overallauc overallauc
a=mean(overallauc)
b=std(overallauc)

