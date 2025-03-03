function checked = check_file_structure(type,file_name)
checked = true;
structure = load(file_name);
switch lower(type)
    case 'cortex'
        if ~isfield(structure, 'Faces') || isempty(structure.Faces)
            checked = false;
        end
        if ~isfield(structure, 'Vertices') || isempty(structure.Vertices)
            checked = false;
        end
        if ~isfield(structure,'VertConn') || isempty(structure.VertConn)
            checked = false;
        end
        if ~isfield(structure, 'VertNormals') || isempty(structure.VertNormals)
            checked = false;
        end
        if ~isfield(structure, 'Curvature') || isempty(structure.Curvature)
            checked = false;
        end
        if ~isfield(structure, 'SulciMap') || isempty(structure.SulciMap)
            checked = false;
        end
        if ~isfield(structure, 'Atlas') || isempty(structure.Atlas) || ~isfield(structure, 'iAtlas') || isempty(structure.iAtlas) || (structure.iAtlas > length(structure.Atlas))
            checked = false;
        elseif ~isfield(structure.Atlas(structure.iAtlas), 'Scouts')
            checked = false;
        end

    case 'leadfield'
        if ~isfield(structure, 'Gain') || isempty(structure.Gain)
            checked = false;
        end
        if ~isfield(structure, 'GridOrient') || isempty(structure.GridOrient)
            checked = false;
        end

    case 'channels'
        if ~isfield(structure, 'Channel') || isempty(structure.Channel)
            checked = false;
        end

    case 'conn_anat'
        if ~isfield(structure, 'Fpt') || isempty(structure.Fpt)
            checked = false;
        end
        if ~isfield(structure, 'parcelIDs') || isempty(structure.parcelIDs)
            checked = false;
        end

    case 'trac'
        if ~isfield(structure, 'parcelIDs') || isempty(structure.parcelIDs)
            checked = false;
        end
        if ~isfield(structure, 'tractLengths') || isempty(structure.tractLengths)
            checked = false;
        end
        if ~isfield(structure, 'tractLengthUnit') || isempty(structure.tractLengthUnit)
            checked = false;
        end

    case 'delay'
        if ~isfield(structure, 'delays') || isempty(structure.delays)
            checked = false;
        end
        if ~isfield(structure, 'delaysMean') || isempty(structure.delaysMean)
            checked = false;
        end
        if ~isfield(structure, 'delaysMedian') || isempty(structure.delaysMedian)
            checked = false;
        end
        if ~isfield(structure, 'delaysMedianSquare') || isempty(structure.delaysMedianSquare)
            checked = false;
        end
        if ~isfield(structure, 'parcelIDs') || isempty(structure.parcelIDs)
            checked = false;
        end
end
end

