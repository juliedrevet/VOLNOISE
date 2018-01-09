%  RUNNER  Runner script for VOLNOISE IRM experiment
%
%  This script runs the VOLNOISE IRM experiment and saves automatically the obtained
%  results in the ../Data folder. This script should be run from its containing
%  folder (which contains all the functions of the VOLNOISE IRM experiment).
%  Original: Valentin Wyart - 09/2015, modified by Julie Drevet 
% =========================================================================

% clear workspace
clear all java
close all hidden
clc

addpath('./Toolboxes/Rand');

% setup random number generator
SetupRand;

% get participant information
prompt={'Subject number (two-digit)','Session number (1/2)','Use CENIR IRM setup? (yes/no)','load raw experiment structure? (enter starting block)'};

argindlg = inputdlg(prompt,'VOLNOISE IRM',1,{'','','no','no'});

if isempty(argindlg)
    error('experiment cancelled!');
end
if isempty(argindlg{1})||isempty(argindlg{2})
    error('experiment cancelled!');
end

hdr         = [];
hdr.subj    = str2num(argindlg{1});
hdr.session = str2num(argindlg{2});
hdr.date    = datestr(now,'yyyymmdd-HHMM');

% use CENIR IRM setup?
%fmri = strcmpi(argindlg{3},'yes'); % strcmpi should work from 2016 on
fmri = strcmp(lower(argindlg{3}),'yes');

% run experiment (or start at specified block)
if ~strcmp(argindlg{4},'no')
    [FileName,PathName] = uigetfile('*.mat','Select subject experiment structure');
    expe_raw = importdata(fullfile(PathName,FileName));
    if isempty(expe_raw)
        error('no experiment structure available');
    end
    start_blck = str2num(argindlg{4});
    [expe,aborted,errmsg] = run_expe(hdr.subj,hdr.session,fmri,expe_raw,start_blck);
    fpath = sprintf('./Data');
    fname = sprintf('VOLNOISE_IRM_S%02d_session%d_%s',hdr.subj,hdr.session,datestr(now,'yyyymmdd-HHMMSS'));
    fname = fullfile(fpath,fname);
    save([fname,'.mat'],'expe');
else
    % run experiment
    [expe,aborted,errmsg] = run_expe(hdr.subj,hdr.session,fmri);
    fpath = sprintf('./Data');
    fname = sprintf('VOLNOISE_IRM_S%02d_session%d_%s',hdr.subj,hdr.session,datestr(now,'yyyymmdd-HHMMSS'));
    fname = fullfile(fpath,fname);
    save([fname,'.mat'],'expe');
end

% add header to experiment structure
expe.hdr = hdr;
expe = orderfields(expe,{'hdr','blck','rslt','clck','stim','logi'});

if ~isempty(errmsg)
    rethrow(errmsg);
end

% save data
datapath = sprintf('./Data/S%02d',hdr.subj);
filename = sprintf('VOLNOISE_IRM_S%02d_session%d_%s',hdr.subj,hdr.session,hdr.date);
if aborted
    filename = [filename,'_aborted'];
end
filename = [filename,'.mat'];
filename = fullfile(datapath,filename);
save(filename,'expe');
