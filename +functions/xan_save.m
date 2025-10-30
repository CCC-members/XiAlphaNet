function [Output] = xan_save(properties,claf,varargin)
import functions.*
import functions.StochasticFISTA.*
import functions.auxx.*
import functions.auxx.BayesOptimization.*
import functions.auxx.CrossSpectrum.*
import functions.auxx.ModelVectorization.*
import tools.*

for i=1:length(varargin)
    eval([inputname(i+2) '= varargin{i};']);
end

switch lower(claf)
    case 'structural'
        structural_path = fullfile(properties.general_params.output_path,'structural');
        if(~isfolder(structural_path))
            mkdir(structural_path);
        end
        save(fullfile(structural_path,'cortex.mat'),'-struct','Cortex');
        save(fullfile(structural_path,'leadfield.mat'),'-struct','Leadfield');
        save(fullfile(structural_path,'parameters.mat'),'-struct','parameters');
        XIALPHANET.Structural.Cortex        = "structural/cortex.mat";
        XIALPHANET.Structural.Leadfield     = "structural/leadfield.mat";
        XIALPHANET.Structural.parameters    = "structural/parameters.mat";
        Output = XIALPHANET;
    case 'create_subject'
        subject_path                            = fullfile(properties.general_params.output_path,SubID);
        if(~isfolder(subject_path))
            mkdir(subject_path);
        end
        Participant.Data = fullfile(strcat(SubID,'_desc-prep_eeg.mat'));
        save(fullfile(subject_path,strcat(SubID,'_desc-prep_eeg.mat')),'-struct','data');
        Output = Participant;
    case 'subject' 
        subject_path                            = fullfile(properties.general_params.output_path,SubID);
        if(~isfolder(subject_path))
            mkdir(subject_path);
        end
        [e,a,~]                                 = x2v(x.Solution);
        freq                                    = parameters.Data.freq;
        % Saving Modulatory Weights of the Anatomical Connectivity and
        % Conduction Delays
        Mod_Weights                             = x.Lambda_DC;
        Participant.Mod_Weights                 = fullfile('Mod_Weights.mat');
        save(fullfile(subject_path,"Mod_Weights.mat"),"Mod_Weights");
        % Saving Parameters of the Xi-Process (Aperiodic Component)
        Xi_estimate.Power                       = e(:,1);    % Xi - Power
        Xi_estimate.Width                       = e(:,2);    % Xi - Width
        Xi_estimate.Exponent                    = e(:,3); % Xi - Exponent
        Participant.Xi_estimate                 = fullfile('Xi_estimate.mat');
        save(fullfile(subject_path,"Xi_estimate.mat"),'-struct',"Xi_estimate");
        % Saving Parameters of the Alpha-Process
        Alpha_estimate.Power                    = a(:,1);    % Alpha - Power
        Alpha_estimate.Width                    = a(:,2);    % Alpha - Width
        Alpha_estimate.Exponent                 = a(:,3); % Alpha - Exponent
        Alpha_estimate.PAF                      = a(:,4);      % Alpha - Peak Frequency
        Participant.Alpha_estimate              = fullfile('Alpha_estimate.mat');
        save(fullfile(subject_path,"Alpha_estimate.mat"),'-struct',"Alpha_estimate");
        % Saving Estimated Condunction Delay Matrix between Neurotracts
        Delay_Matrix                            = parameters.Compact_Model.D;
        Participant.Delay_Matrix                = fullfile('Delay_Matrix.mat');
        save(fullfile(subject_path,"Delay_Matrix.mat"),"Delay_Matrix");
        % Saving Estimated Anatomical Connectivity Matrix between Neurotracts
        Conn_Matrix                             = parameters.Compact_Model.C;
        Participant.Conn_Matrix                 = fullfile('Conn_Matrix.mat');
        save(fullfile(subject_path,"Conn_Matrix.mat"),"Conn_Matrix");
        % %Saving Transfer Function from Source to Scalp
        % Transfer_Function = T;
        % Participant.Transfer_Function = fullfile('Transfer_Function.mat');
        % save(fullfile(subject_path,"Transfer_Function.mat"),"Transfer_Function");
        % Map Solution into Source  Activation and Crosspectrum on the Frequency Domain
        [source_act_cross]                      = eval_source_conn(x.Solution, data.freq,parameters.Model.R,properties,parameters);
        % Saving  Full Activation
        Source_Activations_Full                 = source_act_cross.Activations.Full;
        Participant.Source_Activations_Full     = fullfile('Source_Activations_Full.mat');
        save(fullfile(subject_path,"Source_Activations_Full.mat"),"Source_Activations_Full");
        % Saving Xi-Process Activation
        Source_Activations_Xi                   = source_act_cross.Activations.Xi;
        Participant.Source_Activations_Xi       = fullfile('Source_Activations_Xi.mat');
        save(fullfile(subject_path,"Source_Activations_Xi.mat"),"Source_Activations_Xi");
        % Saving Alpha-Process Activation
        Source_Activations_Alpha                = source_act_cross.Activations.Xi;
        Participant.Source_Activations_Alpha    = fullfile('Source_Activations_Alpha.mat');
        save(fullfile(subject_path,"Source_Activations_Alpha.mat"),"Source_Activations_Alpha");
        % Saving Full Cross-Spectrum in ROIs-Space
        Source_Cross_Full                       = source_act_cross.Cross.Spectra;
        [Source_PSD]                            = 10*log10(Source_Cross_Full);
        Participant.Source_PSD                  = fullfile('Source_PSD.mat');
        save(fullfile(subject_path,"Source_PSD.mat"),"Source_PSD",'freq');

        if(properties.general_params.save_mode.extended) 
            % Saving Full Cross-Spectrum in ROIs-Space
            Source_Cross_Full                   = source_act_cross.Cross.Full;
            Participant.Source_Cross_Full       = fullfile('Source_Cross_Full.mat');
            save(fullfile(subject_path,"Source_Cross_Full.mat"),"Source_Cross_Full");
            % Saving Xi-Process Cross-Spectrum in ROIs-Space
            Source_Cross_Xi                     = source_act_cross.Cross.Xi;
            Participant.Source_Cross_Xi         = fullfile('Source_Cross_Xi.mat');
            save(fullfile(subject_path,"Source_Cross_Xi.mat"),"Source_Cross_Xi");
            % Saving Alpha-Process Cross-Spectrum in ROIs-Space
            Source_Cross_Alpha                  = source_act_cross.Cross.Alpha;
            Participant.Source_Cross_Alpha      = fullfile('Source_Cross_Alpha.mat');
            save(fullfile(subject_path,"Source_Cross_Alpha.mat"),"Source_Cross_Alpha");
        end
        % Saving Status
        Participant.Status                      = "Completed";
        saveJSON(Participant,fullfile(subject_path ,strcat(SubID,'.json')));
        Output = Participant;
    case 'groups'
        group_path = fullfile(properties.general_params.output_path,'groups');
        if(~isfolder(group_path))
            mkdir(group_path);
        end
        disp("-->> Saving Group analysis outputs");
        for i=1:length(properties.model_params.group.age_range)
            age_range = properties.model_params.group.age_range{i};
            age_rage_name = strcat(num2str(age_range{1}),'-',num2str(age_range{2}),'years');
            disp(strcat("---->> Saving group analysis for age range: ",age_rage_name));
            age_range_path = fullfile(group_path,age_rage_name);
            if(~isfolder(age_range_path))
                mkdir(age_range_path);
            end            
            XIALPHANET.Groups.AgeRange(i).Comment = "Some Comment";
            XIALPHANET.Groups.AgeRange(i).Range = age_range;
            disp(strcat("------>> Saving file: ",'file1.mat'));
            save(fullfile(age_range_path,'file1.mat'),'-struct','data1');
            XIALPHANET.Groups.AgeRange(i).FileName = fullfile('groups',age_rage_name,'file1.mat');
            disp(strcat("------>> Saving file: ",'file2.mat'));
            save(fullfile(age_range_path,'file2.mat'),'-struct','data2');            
            XIALPHANET.Groups.AgeRange(i).FileName = fullfile('groups',age_rage_name,'file2.mat');
            disp(strcat("------>> Saving file: ",'file3.mat'));
            save(fullfile(age_range_path,'file3.mat'),'-struct','data3');
            XIALPHANET.Groups.AgeRange(i).FileName = fullfile('groups',age_rage_name,'file3.mat');
        end
        Output = XIALPHANET;        
end
end