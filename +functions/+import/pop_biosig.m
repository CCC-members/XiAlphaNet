% POP_BIOSIG - import data files into EEGLAB using BIOSIG toolbox
%
% Usage:
%   >> OUTEEG = pop_biosig; % pop up window
%   >> OUTEEG = pop_biosig( filename, 'key', val, ...);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'channels'   - [integer array] list of channel indices
%   'blockrange' - [min max] integer range of data blocks to import, in seconds.
%                  Entering [0 3] will import the first three blocks of data.
%                  Default is empty -> import all data blocks. 
%  'importevent' - ['on'|'off'] import events. Default if 'on'.
%  'importannot' - ['on'|'off'] import annotations (EDF+ only). Default if 'on'
%  'importmex'   - ['on'|'off'] import events with Biosig mexSLOAD as an alternative. Default if 'off'
%  'overflow'    - ['on'|'off'] overflow detection. Default is 'off'
%  'uncalibrated' - ['on'|'off'] import uncalibrated data. Default is 'off'
%  'blockepoch'  - ['on'|'off'] force importing continuous data. Default is 'on'
%  'bdfeventmode' - [integer] see bdf2biosig_events function help. Default is 4.
%  'ref'         - [integer] channel index or index(s) for the reference.
%                  Reference channels are not removed from the data,
%                  allowing easy re-referencing. If more than one
%                  channel, data are referenced to the average of the
%                  indexed channels. WARNING! Biosemi Active II data 
%                  are recorded reference-free, but LOSE 40 dB of SNR 
%                  if no reference is used!. If you do not know which
%                  channel to use, pick one and then re-reference after 
%                  the channel locations are read in. {default: none}.
%                  For more information see http://www.biosemi.com/faq/cms&drl.htm
%  'refoptions'  - [Cell] Option for the pop_reref function. Default is to 
%                  remove the reference channel if there is one of them and to 
%                  keep it if there are several of them from the graphic
%                  interface. From the command line default option is to 
%                  keep the reference channel.
%  'rmeventchan' - ['on'|'off'] remove event channel after event 
%                  extraction. Default is 'on'.
%  'memorymapped' - ['on'|'off'] import memory mapped file (useful if 
%                  encountering memory errors). Default is 'off'.
%
% Outputs:
%   OUTEEG   - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Oct. 29, 2003-
%
% Note: BIOSIG toolbox must be installed. Download BIOSIG at 
%       http://biosig.sourceforge.net
%       Contact a.schloegl@ieee.org for troubleshooting using BIOSIG.

% Copyright (C) 2003 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [ALLEEG, dat] = pop_biosig(filename, varargin)
import functions.*
import functions.import.*
ALLEEG = [];
options = varargin;
saveData = false;
  
% decode input parameters
% -----------------------
g = finputcheck( options, { 'blockrange'   'integer' [0 Inf]    [];
                            'channels'     'integer' [0 Inf]    [];
                            'ref'          'integer' [0 Inf]    [];
                            'bdfeventmode' 'integer' [0 Inf]    4;
                            'refoptions'   'cell'    {}             { 'keepref' 'on' };
                            'rmeventchan'  'string'  { 'on';'off' } 'on';
                            'overflow'     'string'  { 'on';'off' } 'off';
                            'uncalibrated' 'string'  { 'on';'off' } 'off';
                            'importevent'  'string'  { 'on';'off' } 'on';
                            'importannot'  'string'  { 'on';'off' } 'on';
                            'importmex'   'string'  { 'on';'off' }  'off';
                            'memorymapped' 'string'  { 'on';'off' } 'off';
                            'blockepoch'   'string'  { 'on';'off' } 'off' }, 'pop_biosig');
if ischar(g), error(g); end
if ~iscell(filename)
    filename = { filename };
end

for iFile = 1:length(filename)
    % import data
    % -----------
    if ~exist(filename{iFile})
        error('File not found %s', filename{iFile})
    end
    EEG = eeg_emptyset;
    [dat, DAT, interval] = readfile(filename{iFile}, [], g.blockrange, g.memorymapped, g.bdfeventmode, g.overflow, g.uncalibrated);
    
    if strcmpi(g.blockepoch, 'off')
        dat.NRec = 1;
    end
    EEG = biosig2eeglab(dat, DAT, interval, g.channels, strcmpi(g.importevent, 'on'), strcmpi(g.importannot, 'on'));
    
    if strcmpi(g.rmeventchan, 'on') && strcmpi(dat.TYPE, 'BDF') && isfield(dat, 'BDF')
        if size(EEG.data,1) >= dat.BDF.Status.Channel
            disp('Removing event channel...');
            EEG.data(dat.BDF.Status.Channel,:) = [];
            if ~isempty(EEG.chanlocs) && length(EEG.chanlocs) >= dat.BDF.Status.Channel
                EEG.chanlocs(dat.BDF.Status.Channel) = [];
            end
        end
        EEG.nbchan = size(EEG.data,1);
    end
    
    % rerefencing
    % -----------
    if ~isempty(g.ref)
        disp('Re-referencing...');
        EEG = pop_reref(EEG, g.ref, g.refoptions{:});
    end
    
    % test if annotation channel is present
    % -------------------------------------
    if isfield(dat, 'EDFplus') && strcmpi(g.importannot, 'on')
        tmpfields = fieldnames(dat.EDFplus);
        for ind = 1:length(tmpfields)
            tmpdat = getfield(dat.EDFplus, tmpfields{ind});
            if length(tmpdat) == EEG.pnts
                EEG.data(end+1,:) = tmpdat;
                EEG.nbchan        = EEG.nbchan+1;
                if ~isempty(EEG.chanlocs)
                    EEG.chanlocs(end+1).labels = tmpfields{ind};
                end
            end
        end
    end
    
    % import using Biosig mexSLOAD method (Cedric edits 2/23/2021)
    if strcmpi(g.importmex, 'on')
        [s,HDR] = mexSLOAD(filename);
        
        %Get correct event names contained in CodeDesc
        num_ev_type = unique(HDR.EVENT.TYP);
        num_ev_name = unique(HDR.EVENT.CodeDesc);
        if ~isempty(HDR.EVENT.CodeDesc) && length(num_ev_type) == length(num_ev_name)
            for iEvent = 1:length(HDR.EVENT.TYP)
                EEG.event(iEvent).type = char(HDR.EVENT.CodeDesc(HDR.EVENT.TYP(iEvent)));
                EEG.event(iEvent).latency = HDR.EVENT.POS(iEvent);
            end
        else
            for iEvent = 1:length(HDR.EVENT.TYP)
                EEG.event(iEvent).type = HDR.EVENT.TYP(iEvent);
                EEG.event(iEvent).latency = HDR.EVENT.POS(iEvent);
                EEG.event(iEvent).urevent = iEvent;
            end
            warning('Inconsistency between event types and names or event names were not found.');
            warning('Check your events or convert your data format with EDFBrowser: https://www.teuniz.net/edfbrowser/');
        end
    end
    
    % check and store data
    % --------------------
%     EEG = eeg_checkset(EEG, 'makeur');   % Make EEG.urevent field
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = eeg_checkset(EEG);
    if saveData
        EEG = pop_saveset(EEG, [ filename{iFile}(1:end-4) '.set' ]);
        EEG(1).saved = 'justloaded';
    end
    if iFile == 1
        ALLEEG = EEG;
    else
        ALLEEG(iFile) = EEG;
    end
    
end

% history
% -------
if length(filename) == 1
    command = sprintf('EEG = pop_biosig(''%s''', filename{1});
else
    command = sprintf('EEG = pop_biosig(%s', vararg2str(filename));
end
if isempty(options)
    command = [ command ');' ];
else
    command = [ command sprintf(', %s);', vararg2str(options)) ];
end

% Checking if str2double is on top of the path
biosigpathlast;

% ---------
% read data
% ---------
function [dat, DAT, interval] = readfile(filename, channels, blockrange, memmapdata, bdfeventmode, overflow, uncalibrated)

if isempty(channels), channels = 0; end
strmode = '';
if strcmpi(overflow, 'off')
    strmode = 'OVERFLOWDETECTION:OFF';
end
if strcmpi(uncalibrated, 'on')
    if ~isempty(strmode)
        strmode = [strmode ';'];
    end
    strmode = [ strmode 'UCAL'];
end
if ~isequal(bdfeventmode, 4)
    if ~isempty(strmode)
        strmode = [strmode ';'];
    end
    strmode = [ strmode 'BDF:[' num2str(bdfeventmode) ']' ];
end
fprintf('sopen mode is "%s"\n', strmode);
dat = sopen(filename, 'r', channels, strmode);

if strcmpi(memmapdata, 'off')
    fprintf('Reading data in %s format...\n', dat.TYPE);

    if ~isempty(blockrange)
        newblockrange    = blockrange;
%         newblockrange    = newblockrange*dat.Dur;    
        DAT=sread(dat, newblockrange(2)-newblockrange(1), newblockrange(1));
    else 
        DAT=sread(dat, Inf);% this isn't transposed in original!!!!!!!!
        newblockrange    = [];
    end
    sclose(dat);
else
    fprintf('Reading data in %s format (file will be mapped to memory so this may take a while)...\n', dat.TYPE);
    inc = ceil(250000/(dat.NS*dat.SPR)); % 1Mb block
    
    if isempty(blockrange), blockrange = [0 dat.NRec]; end
    blockrange(2) = min(blockrange(2), dat.NRec);
    allblocks = [blockrange(1):inc:blockrange(end)];
    count = 1;
    for bind = 1:length(allblocks)-1
        TMPDAT=sread(dat, (allblocks(bind+1)-allblocks(bind))*dat.Dur, allblocks(bind)*dat.Dur);
        if bind == 1
            DAT = mmo([], [size(TMPDAT,2) (allblocks(end)-allblocks(1))*dat.SPR]);
        end
        DAT(:,count:count+length(TMPDAT)-1) = TMPDAT';
        count = count+length(TMPDAT);
    end
    sclose(dat);
end

if ~isempty(blockrange)
     interval(1) = blockrange(1) * dat.SampleRate(1) + 1;
     interval(2) = blockrange(2) * dat.SampleRate(1);
else interval = [];
end
