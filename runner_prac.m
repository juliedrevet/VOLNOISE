% clear workspace
clear all java
close all hidden
clc

addpath ./Toolboxes/Rand

% setup random number generator
SetupRand;

prompt={'Subject number (two-digit)'};

argindlg = inputdlg(prompt,'VOLNOISE IRM',1,{''});

if isempty(argindlg)
    error('experiment cancelled!');
end
if isempty(argindlg{1})
    error('experiment cancelled!');
end

hdr      = [];
hdr.subj = str2num(argindlg{1});
hdr.date = datestr(now,'yyyymmdd-HHMM');

% generation of practice blocks
nprac = 2; % number of practice blocks

expe_prac = gen_prac(nprac);

% add header to experiment structure
expe_prac.hdr = hdr;
expe_prac     = orderfields(expe_prac,{'hdr','blck'});

[expe_prac,aborted] = run_prac(hdr.subj,expe_prac,true);

% save data
foldname = sprintf('./Data/S%02d',hdr.subj);
filename = sprintf('VOLNOISE_IRM_S%02d_training_%s',hdr.subj,hdr.date);
if aborted
    filename = [filename,'_aborted'];
end
filename = [filename,'.mat'];
filename = fullfile(foldname,filename);
save(filename,'expe_prac');

