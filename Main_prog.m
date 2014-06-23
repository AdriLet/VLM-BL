clear all


tic
 
%%  Load geometry
load('Géométrie/HeadF2M.mat');
load('Géométrie/MainF2M.mat');


%% Calculation of potential flow
i=1;
for awa=10:30
clearvars -except awa HeadF2M MainF2M results i



result=run_vlm(4,5,2,awa,HeadF2M,MainF2M);
Head=result(1);
Main=result(2);

%% Calculation of Streamlines

result=StreamLines([Head Main]);
Head=result(1);
Main=result(2);


%% Calculation of Boundary layer


Head=SailBoundaryLayer(Head);
Main=SailBoundaryLayer(Main);



PlotSails([Head Main]);

results.Cl(i)=Main.outputs.Cl;
results.AWA(i)=Main.AWA;
results.Cdi(i)=Main.outputs.Cdi;

i=i+1;

end;


toc
