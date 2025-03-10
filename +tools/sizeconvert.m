function[sizes_c]=sizeconvert(sizes,rounding)
%SIZECONVERT Convert size from bytes to larger units of size as relevant.
%
%  SIZECONVERT(SIZES,ROUNDING) converts a numerical matrix containing sizes
%  in bytes to larger units of size as relevant. The output is a cell with
%  the first column containing the converted sizes and the second column
%  containing the abbreviated converted units.
%
%  By default, the converted sizes are rounded to the nearest integer for
%  readability, although this can be disabled.
%
%  OPTIONAL ARGUMENTS:
%    ROUNDING:
%      1: Round converted sizes to nearest integer (default).
%      0: Do not round converted sizes to nearest integer.
%
%  EXAMPLES:
%    sizeconvert([534 1019 1047 16854 1724416 189340868 403296067584])
%      returns {534,'B';1019,'B';1,'KB';16,'KB';2,'MB';181,'MB';376,'GB'}
%    sizeconvert([437,5457, 4534405],0)
%      returns {437,'B';5.3291,'KB';4.3243,'MB'}
%
%  VERSION DATE: 2005.04.24
%  MATLAB VERSION: 7.0.1.24704 (R14) Service Pack 1
%
%{
REVISION HISTORY:
2005.04.24: Made minor updates to comments.
2005.04.12: Original release.
KEYWORDS:
size, size conversion
byte, kilobyte, megabyte, gigabyte, terabyte, petabyte
B, KB, MB, GB, TB, PB
%}
%--------------------------------------------------------------------------
%% Set value for ROUNDING if missing.
if nargin==1
    rounding=1;
end
%--------------------------------------------------------------------------
%% Convert sizes.
for count=1:numel(sizes)
    
    size_cur=sizes(count);
    
    if size_cur<1024^1
        sizes_c_cur=size_cur/(1024^0);
        sizes_c_unit_cur='B';
    elseif size_cur<1024^2
        sizes_c_cur=size_cur/(1024^1);
        sizes_c_unit_cur='KB';
    elseif size_cur<1024^3
        sizes_c_cur=size_cur/(1024^2);
        sizes_c_unit_cur='MB';
    elseif size_cur<1024^4
        sizes_c_cur=size_cur/(1024^3);
        sizes_c_unit_cur='GB';
    elseif size_cur<1024^5
        sizes_c_cur=size_cur/(1024^4);
        sizes_c_unit_cur='TB';
    elseif size_cur<1024^6
        sizes_c_cur=size_cur/(1024^5);
        sizes_c_unit_cur='PB';
    else
        error('Unsupported size')
    end
    
    sizes_c(count,1)={sizes_c_cur};
    sizes_c(count,2)={sizes_c_unit_cur};
    
end
%--------------------------------------------------------------------------
%% Round converted sizes if ROUNDING is enabled.
if rounding==1
    sizes_c_sizes=sizes_c(:,1);
    sizes_c_sizes=cell2mat(sizes_c_sizes);
    sizes_c_sizes=round(sizes_c_sizes);
    sizes_c(:,1)=num2cell(sizes_c_sizes);
end
%--------------------------------------------------------------------------

