% read pressure in input file

file = './yieldStress/postProcessing/inletPressure/111601.7/surfaceFieldValue.dat';

opts = detectImportOptions(file);
data = readmatrix(file, opts);

t = data(:,1);
p = data(:,2);



figure();
plot(t/3600, p, 'b.')
xlabel('time [h]');
ylabel('pressure [Pa]')