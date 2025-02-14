function [Participant] = xan_save(properties,SubID,claf,varargin)

for i=1:length(varargin)
    eval([inputname(i+3) '= varargin{i};']);
end

switch lower(claf)
    case 'create_subject'
        subject_path                    = fullfile(properties.general_params.output_path,SubID);
        if(~isfolder(subject_path))
            mkdir(subject_path);
        end
        Participant.Data = fullfile(strcat(SubID,'_desc-prep_eeg.mat'));
        save(fullfile(subject_path,strcat(SubID,'_desc-prep_eeg.mat')),'-struct','data');

    case 'subject' 
        subject_path                    = fullfile(properties.general_params.output_path,SubID);
        if(~isfolder(subject_path))
            mkdir(subject_path);
        end
        [e,a,~] = x2v(x.Solution);
        % Saving Modulatory Weights of the Anatomical Connectivity and
        % Conduction Delays
        Mod_Weights = x.Lambda_DC;
        Participant.Mod_Weights = fullfile('Mod_Weights.mat');
        save(fullfile(subject_path,"Mod_Weights.mat"),"Mod_Weights");
        % Saving Parameters of the Xi-Process (Aperiodic Component)
        Xi_estimate.Power = e(:,1);    % Xi - Power
        Xi_estimate.Width = e(:,2);    % Xi - Width
        Xi_estimate.Exponent = e(:,3); % Xi - Exponent
        Participant.Xi_estimate = fullfile('Xi_estimate.mat');
        save(fullfile(subject_path,"Xi_estimate.mat"),'-struct',"Xi_estimate");
        % Saving Parameters of the Alpha-Process
        Alpha_estimate.Power = a(:,1);    % Alpha - Power
        Alpha_estimate.Width = a(:,2);    % Alpha - Width
        Alpha_estimate.Exponent = a(:,3); % Alpha - Exponent
        Alpha_estimate.PAF = a(:,4);      % Alpha - Peak Frequency
        Participant.Alpha_estimate = fullfile('Alpha_estimate.mat');
        save(fullfile(subject_path,"Alpha_estimate.mat"),'-struct',"Alpha_estimate");
        % Saving Estimated Condunction Delay Matrix between Neurotracts
        Delay_Matrix = x.Lambda_DC(1)*parameters.Compact_Model.D;
        Participant.Delay_Matrix = fullfile('Delay_Matrix.mat');
        save(fullfile(subject_path,"Delay_Matrix.mat"),"Delay_Matrix");
        % Saving Estimated Anatomical Connectivity Matrix between Neurotracts
        Conn_Matrix = x.Lambda_DC(2)*parameters.Compact_Model.C;
        Participant.Conn_Matrix = fullfile('Conn_Matrix.mat');
        save(fullfile(subject_path,"Conn_Matrix.mat"),"Conn_Matrix");
        %Saving Transfer Function from Source to Scalp
        Transfer_Function = T;
        Participant.Transfer_Function = fullfile('Transfer_Function.mat');
        save(fullfile(subject_path,"Transfer_Function.mat"),"Transfer_Function");

        if(properties.general_params.save_mode.extended)
            % Map Solution into Source  Activation and Crosspectrum on the Frequency Domain
            [source_act_cross] = eval_source_conn(x.Solution, data.freq,T,parameters.Model.K,parameters.Model.R,properties);
            % Saving  Full Activation
            Source_Activations_Full = source_act_cross.Activations.Full;
            Participant.Source_Activations_Full = fullfile('Source_Activations_Full.mat');
            save(fullfile(subject_path,"Source_Activations_Full.mat"),"Source_Activations_Full");
            % Saving Xi-Process Activation
            Source_Activations_Xi = source_act_cross.Activations.Xi;
            Participant.Source_Activations_Xi = fullfile('Source_Activations_Xi.mat');
            save(fullfile(subject_path,"Source_Activations_Xi.mat"),"Source_Activations_Xi");
            % Saving Alpha-Process Activation
            Source_Activations_Alpha = source_act_cross.Activations.Xi;
            Participant.Source_Activations_Alpha = fullfile('Source_Activations_Alpha.mat');
            save(fullfile(subject_path,"Source_Activations_Alpha.mat"),"Source_Activations_Alpha");
            % Saving Full Cross-Spectrum in ROIs-Space
            Source_Cross_Full = source_act_cross.Cross.Full;
            Participant.Source_Cross_Full = fullfile('Source_Cross_Full.mat');
            save(fullfile(subject_path,"Source_Cross_Full.mat"),"Source_Cross_Full");
            % Saving Xi-Process Cross-Spectrum in ROIs-Space
            Source_Cross_Xi = source_act_cross.Cross.Xi;
            Participant.Source_Cross_Xi = fullfile('Source_Cross_Xi.mat');
            save(fullfile(subject_path,"Source_Cross_Xi.mat"),"Source_Cross_Xi");
            % Saving Alpha-Process Cross-Spectrum in ROIs-Space
            Source_Cross_Alpha = source_act_cross.Cross.Alpha;
            Participant.Source_Cross_Alpha = fullfile('Source_Cross_Alpha.mat');
            save(fullfile(subject_path,"Source_Cross_Alpha.mat"),"Source_Cross_Alpha");
        end
        % Saving Status
        Participant.Status = "Completed";
        saveJSON(Participant,fullfile(subject_path ,strcat(SubID,'.json')));
end
end