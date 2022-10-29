close('all')


Directory = '../../../data/LargeScale/B13289O14-DH1-01463/Day1-09_10-12-21/DatData/Clipped/B13289O14-DH1-Rec9inter';

Par =  LoadXml([Directory 'AC.xml']);

d = dir(strcat(Directory, '*.lfp'));%DC-LP30Hz-Notch50-100Hz.dat');

LFPFs = 651.04166667;
DownFact = 100;